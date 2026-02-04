"""
================================================================================
File:        cov-est.smk
Description: Mean coverage and relative abundance estimation of vOTUs
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/coverage-estimation
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    snakemake -s cov-est.smk --cores 48

================================================================================
"""

configfile: "config.yaml"
REF = config["reference_fasta"]

import pandas as pd

samples_df = pd.read_csv("samples.tsv", sep="\t").set_index("sample_id", drop=False)

def get_R1(wildcards):
    return samples_df.loc[wildcards.sample, "R1_path"]
def get_R2(wildcards):
    return samples_df.loc[wildcards.sample, "R2_path"]

rule all:
    input:
        "results/" + config["output_prefix"] + "_abundance.tsv"

# --- Map reads to reference ---

rule bwa_index:
    input:
        REF
    output:
        multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        ref = REF,
        idx = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        r1 = get_R1,
        r2 = get_R2
    output:
        temp("results/{sample}.sam")
    threads: 8
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} > {output}"

# --- Filter reads ---

rule filter_and_sort:
    input:
        "results/{sample}.sam"
    output:
        temp("results/{sample}.sorted.bam")
    threads: 8
    shell:
        "samtools view -h -F 256 -u {input} | samtools sort -@ {threads} -o {output}"

rule index_bam:
    input:
        "results/{sample}.sorted.bam"
    output:
        temp("results/{sample}.sorted.bam.bai")
    threads: 8
    shell:
        "samtools index -@ {threads} {input} {output}"

rule filter_bam:
    input:
        bam = "results/{sample}.sorted.bam",
        bai = "results/{sample}.sorted.bam.bai"
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

# --- Generate coverage and depth statistics ---

rule coverage:
    input:
        "results/{sample}_flt.bam"
    output:
        "results/{sample}_coverage.tsv"
    shell:
        "samtools coverage {input} > {output}"

rule depth:
    input:
        "results/{sample}_flt.bam"
    output:
        "results/{sample}_depth.tsv"
    shell:
        "samtools depth {input} > {output}"

# --- Compute clipped mean coverage ---

rule get_clipped_coverage:
    input:
        cov = "results/{sample}_coverage.tsv",
        depth = "results/{sample}_depth.tsv"
    output:
        "results/{sample}_clippedcov.tsv"
    params:
        prefix = config["ref_header_prefix"]
    shell:
        "python ../scripts/get_clipped_coverage.py -c {input.cov} -d {input.depth} -s {wildcards.sample} -p {params.prefix} -o {output}"

# --- Detection and abundance estimation ---

rule detection_abundance:
    input:
        "results/{sample}_clippedcov.tsv"
    output:
        "results/{sample}_abundance.tsv"
    params:
        covthresh = config["covthresh"]
    run:
        df = pd.read_csv(input[0], sep='\t')

        df = df[df.covbreadth >= params.covthresh]

        df['sum_covmean'] = df.groupby('sample_id').covmean.transform('sum')
        df['perc_abun'] = df.covmean/df['sum_covmean']*100
        df = df.drop(columns='sum_covmean')

        df.to_csv(output[0], sep='\t', index=False)

# --- Aggregate results ---

rule combine_results:
    input:
        expand("results/{sample}_abundance.tsv", sample=samples_df.index)
    output:
        "results/" + config["output_prefix"] + "_abundance.tsv"
    run:
        dfs = [pd.read_csv(tsv, sep='\t') for tsv in input]
        df = pd.concat(dfs)
        df.to_csv(output[0], sep='\t', index=False)
