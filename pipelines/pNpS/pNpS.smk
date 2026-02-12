"""
================================================================================
File:        pNpS.smk
Description: Annotation of possible synonymous and nonsynonymous SNPs in genes;
             read mapping to a set of input vOTUs and decoy MAGs followed by
             BAM filtering, variant calling, and pN/pS analysis.
Authors:     Windsor Koh, Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/pNpS
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    snakemake -s pNpS.smk --cores 48

================================================================================
"""

configfile: "config.yaml"
REF_FASTA = config["reference_fasta"]
DECOY_FASTA = config["decoy_fasta"]
OUT_PREFIX = config["output_prefix"]
REF_PLUS_DECOY = "results/ref_plus_decoy.fna"

import pandas as pd
import os
import sys

REF_PREFIX = os.path.splitext(os.path.basename(REF_FASTA))[0]

samples_df = pd.read_csv("samples.tsv", sep="\t").set_index("sample_id", drop=False)

def get_R1(wildcards):
    return samples_df.loc[wildcards.sample, "R1_path"]
def get_R2(wildcards):
    return samples_df.loc[wildcards.sample, "R2_path"]

rule all:
    input:
        "results/" + OUT_PREFIX + "_vcf.tsv",
        "results/" + OUT_PREFIX + "_gene_pNpS.tsv"

# --- Gene identification and synonymous/nonsynonymous site counts ---

rule genomad:
    input: REF_FASTA
    output:
        dir = directory(f"results/{REF_PREFIX}_ganno"),
        gene_file = f"results/{REF_PREFIX}_ganno/{REF_PREFIX}_annotate/{REF_PREFIX}_genes.tsv"
    params:
        db = config["genomad_db"]
    threads: 16
    shell:
        "genomad annotate -t {threads} {input} {output.dir} {params.db}"

rule possible_site_counts:
    input:
        gene_file = f"results/{REF_PREFIX}_ganno/{REF_PREFIX}_annotate/{REF_PREFIX}_genes.tsv",
        fasta_file = REF_FASTA
    output:
        processed_gene_file = f"results/{REF_PREFIX}_genes_processed.tsv"
    run:
        sys.path.append(os.path.abspath(os.path.join(workflow.basedir, "..")))
        from scripts.pNpS import preprocess_gene_file_possible_sites
        df = preprocess_gene_file_possible_sites(input.gene_file, input.fasta_file)
        df.to_csv(output.processed_gene_file, sep='\t', index=False)

# --- Map reads to reference and decoys ---

rule combine_fasta:
    input:
        ref = REF_FASTA,
        decoy = DECOY_FASTA
    output:
        temp(REF_PLUS_DECOY)
    shell:
        "cat {input.ref} {input.decoy} > {output}"

rule bwa_index:
    input:
        REF_PLUS_DECOY
    output:
        temp(multiext(REF_PLUS_DECOY, ".amb", ".ann", ".bwt", ".pac", ".sa"))
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        ref = REF_PLUS_DECOY,
        idx = multiext(REF_PLUS_DECOY, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        r1 = get_R1,
        r2 = get_R2
    output:
        temp("results/{sample}_ref_plus_decoy.sam")
    threads: 8
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} > {output}"

# --- Filter reads ---

rule get_ref_bed:
    input:
        REF_FASTA
    output:
        fai = REF_FASTA + ".fai",
        bed = "results/ref.bed"
    shell:
        """
        samtools faidx {input}
        awk 'BEGIN{{OFS="\\t"}} {{print $1, 0, $2}}' {output.fai} > {output.bed}
        """

rule filter_and_sort:
    input:
        sam = "results/{sample}_ref_plus_decoy.sam",
        bed = "results/ref.bed"
    output:
        temp("results/{sample}.sorted.bam")
    threads: 8
    shell:
        """
        samtools view -h -L {input.bed} -F 256 -u {input.sam} | \
        samtools sort -@ {threads} -o {output}
        """

rule filter_bam:
    input:
        bam = "results/{sample}.sorted.bam"
    output:
        "results/{sample}_flt.bam"
    params:
        min_id = config["min_id"],
        min_cov = config["min_cov"],
        check_proper_pair = config["check_proper_pair"]
    shell:
        """
        python ../scripts/bamfilter.py \
            -i {input.bam} -o {output} \
            --min_id {params.min_id} --min_cov {params.min_cov} \
            --check_proper_pair {params.check_proper_pair}
        """

# --- Generate coverage statistics ---

rule coverage:
    input:
        "results/{sample}_flt.bam"
    output:
        "results/{sample}_coverage.tsv"
    shell:
        "samtools coverage {input} > {output}"

# --- Call variants using LoFreq ---

rule index_flt_bam:
    input:
        "results/{sample}_flt.bam"
    output:
        temp("results/{sample}_flt.bam.bai")
    threads: 2
    shell:
        "samtools index -@ {threads} {input} {output}"

rule lofreq:
    input:
        bam = "results/{sample}_flt.bam",
        bai = "results/{sample}_flt.bam.bai",
        ref = REF_FASTA
    output:
        "results/{sample}.vcf"
    threads: 8
    shell:
        "lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output} {input.bam}"

# --- Compute pN/pS statistics ---

rule run_pNpS:
    input:
        cov_file = "results/{sample}_coverage.tsv",
        vcf_file = "results/{sample}.vcf",
        processed_gene_file = f"results/{REF_PREFIX}_genes_processed.tsv",
        fasta_file = REF_FASTA
    output:
        vcf_out = temp("results/{sample}_vcf.tsv"),
        gene_pNpS = temp("results/{sample}_gene_pNpS.tsv")
    params:
        covthresh = config["covthresh"]
    shell:
        """
        python ../scripts/pNpS.py --sample {wildcards.sample} --gene_file {input.processed_gene_file} \
        --fasta_file {input.fasta_file} --cov_file {input.cov_file} --vcf_file {input.vcf_file} \
        --cov_threshold {params.covthresh} --vcf_out {output.vcf_out} --gene_pNpS {output.gene_pNpS}
        """

# --- Aggregate results ---

def concat_and_save(input_files, output_file):
    dfs = [pd.read_csv(f, sep='\t') for f in input_files]
    df = pd.concat(dfs)
    df.to_csv(output_file, sep='\t', index=False)

rule merge_results:
    input:
        vcf_in = expand("results/{sample}_vcf.tsv", sample=samples_df.index),
        gene_in = expand("results/{sample}_gene_pNpS.tsv", sample=samples_df.index)
    output:
        vcf_out = "results/" + OUT_PREFIX + "_vcf.tsv",
        gene_out = "results/" + OUT_PREFIX + "_gene_pNpS.tsv"
    run:
        concat_and_save(input.vcf_in, output.vcf_out)
        concat_and_save(input.gene_in, output.gene_out)
