"""
================================================================================
File:        hic.smk
Description: Maps Hi-C reads to metagenomic assemblies, filters and name-sorts
             BAM files, computes number of intra and inter-MAG linkages to infer
             noise model parameters, extracts number of virus-MAG linkages,
             computes adjusted p-values.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/HiC-host-association
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    snakemake -s hic.smk --cores 48

================================================================================
"""

configfile: "config.yaml"

import pandas as pd

samples_df = pd.read_csv("samples.tsv", sep="\t").set_index("sample_id", drop=False)

SAMPLES = samples_df.index
def get_contigs(wildcards):
    return samples_df.loc[wildcards.sample, "contigs_path"]
def get_HiC_R1(wildcards):
    return samples_df.loc[wildcards.sample, "HiC_R1_path"]
def get_HiC_R2(wildcards):
    return samples_df.loc[wildcards.sample, "HiC_R2_path"]
def get_mag_contig_info(wildcards):
    return samples_df.loc[wildcards.sample, "mag_contig_list"]
def get_viral_contig_info(wildcards):
    return samples_df.loc[wildcards.sample, "viral_contig_list"]

import numpy as np
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

def in_region(pos_start, pos_end, region_start, region_end):
    return 1 if pos_start >= region_start and pos_end <= region_end else 0

def binomial_pval(N, p, observed):
    return binom.sf(observed - 1, N, p)

rule all:
    input:
        "results/virus-mag_HiC.tsv"

# --- Preprocess Hi-C reads ---

rule fastp:
    input:
        r1 = get_HiC_R1,
        r2 = get_HiC_R2
    output:
        r1_trimmed = "results/{sample}/{sample}_HiC.1.fastq.gz",
        r2_trimmed = "results/{sample}/{sample}_HiC.2.fastq.gz",
        json = "results/{sample}/{sample}_HiC.json",
        html = "results/{sample}/{sample}_HiC.html"
    threads: 8
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_trimmed} -O {output.r2_trimmed} \
              -j {output.json} -h {output.html} -w {threads}
        """

# --- Map Hi-C reads to metagenomic assemblies ---

rule bwa_index:
    input:
        get_contigs
    output:
        ref = "results/{sample}/{sample}.fna",
        idx = multiext("results/{sample}/{sample}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell:
        """
        ln -rs {input} {output.ref}
        bwa index {output.ref}
        """

rule bwa_mem:
    input:
        ref = "results/{sample}/{sample}.fna",
        idx = multiext("results/{sample}/{sample}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        r1_trimmed = "results/{sample}/{sample}_HiC.1.fastq.gz",
        r2_trimmed = "results/{sample}/{sample}_HiC.2.fastq.gz"
    output:
        temp("results/{sample}/{sample}_HiC.sam")
    threads: 8
    shell:
        "bwa mem -5SP -t {threads} {input.ref} {input.r1_trimmed} {input.r2_trimmed} | samblaster > {output}"

# --- Filter reads ---

rule view_sort:
    input:
        "results/{sample}/{sample}_HiC.sam"
    output:
        temp("results/{sample}/{sample}_HiC.sorted.bam")
    threads: 8
    shell:
        "samtools view -h -F 3340 -u {input} | samtools sort -n -@ {threads} -o {output}"

rule filter_bam:
    input:
        bam = "results/{sample}/{sample}_HiC.sorted.bam"
    output:
        "results/{sample}/{sample}_HiC_flt.bam"
    params:
        min_id = config["min_id"],
        min_cov = config["min_cov"]
    shell:
        """
        python ../scripts/bamfilter.py \
            -i {input.bam} -o {output} \
            --min_id {params.min_id} --min_cov {params.min_cov} \
            --check_proper_pair 0
        """

# --- Construct Hi-C noise model for MAGs ---

rule compile_linkages_between_contigs:
    input:
        "results/{sample}/{sample}_HiC_flt.bam"
    output:
        "results/{sample}/{sample}_numlinkages.tsv"
    shell:
        "python ../scripts/bam_to_numlinkages.py -i {input} -o {output}"

rule compile_linkages_between_mags:
    input:
        linkages = "results/{sample}/{sample}_numlinkages.tsv",
        magcontigs = get_mag_contig_info
    output:
        "results/{sample}/{sample}_maglinkages.tsv"
    run:
        df = pd.read_csv(input.linkages, sep='\t')
        df_c = pd.read_csv(input.magcontigs, sep='\t')

        df = df.merge(df_c.add_suffix('1'), on='contig_id1').merge(df_c.add_suffix('2'), on='contig_id2')

        m1 = df['mag_id1']
        m2 = df['mag_id2']
        df['mag_A'] = m1.where(m1 < m2, m2)
        df['mag_B'] = m2.where(m1 < m2, m1)
        df1 = df.groupby(['mag_A', 'mag_B'], as_index=False)['num_linkages'].sum()

        df1 = df1.rename(columns={'mag_A': 'mag_id1', 'mag_B': 'mag_id2'}).sort_values('num_linkages', ascending=False)

        df1.to_csv(output[0], sep='\t', index=False)

rule solve_HiC_model:
    input:
        "results/{sample}/{sample}_maglinkages.tsv"
    output:
        "results/{sample}/{sample}_mag_HiC_stats.tsv"
    params:
        e_init = config["hic_model"]["e_init"],
        max_iter = config["hic_model"]["max_iter"],
        tol = config["hic_model"]["tol"]
    shell:
        """
        python ../scripts/solve_hic_model.py \
            -i {input} -o {output} \
            --e_init {params.e_init} --max_iter {params.max_iter} --tol {params.tol}
        """

# --- Compute virus-MAG linkages ---

rule get_v_and_binlist:
    input:
        vcontigs = get_viral_contig_info,
        magcontigs = get_mag_contig_info
    output:
        vlist = "results/{sample}/{sample}_viral.list",
        v_and_maglist = "results/{sample}/{sample}_viral_and_mags.list"
    run:
        v = set(pd.read_csv(input.vcontigs, sep='\t', usecols=['contig_id'])['contig_id'])
        with open(output.vlist, "w") as f:
            f.write("\n".join(sorted(v)) + "\n")
        m = set(pd.read_csv(input.magcontigs, sep='\t', usecols=['contig_id'])['contig_id'])
        combined = v | m
        with open(output.v_and_maglist, "w") as f:
            f.write("\n".join(sorted(combined)) + "\n")

rule compile_linkage_positions_between_viral_and_mags:
    input:
        bam = "results/{sample}/{sample}_HiC_flt.bam",
        vlist = "results/{sample}/{sample}_viral.list",        
        v_and_maglist = "results/{sample}/{sample}_viral_and_mags.list"
    output:
        linkage_positions = "results/{sample}/{sample}_linkage_positions_of_viral_mag.tsv"
    shell:
        "python ../scripts/bam_to_mapping_positions.py -i {input.bam} -l1 {input.vlist} -l2 {input.v_and_maglist} -o {output.linkage_positions}"

rule annotate_viral_regions:
    input:
        linkage_positions = "results/{sample}/{sample}_linkage_positions_of_viral_mag.tsv",
        vcontigs = get_viral_contig_info,
        magcontigs = get_mag_contig_info
    output:
        linkage_positions_annotated = "results/{sample}/{sample}_linkage_positions_of_viral_mag_annot.tsv"
    run:
        df = pd.read_csv(input.linkage_positions, sep='\t')
        df_v = pd.read_csv(input.vcontigs, sep='\t').add_suffix('1')
        df_c = pd.read_csv(input.magcontigs, sep='\t').add_suffix('2')

        df = df_v.merge(df, on='contig_id1', how='left').merge(df_c, on='contig_id2', how='left')

        df["in_viral_region1"] = df.apply(lambda row: in_region(row.pos1_start, row.pos1_end, row.start_bp1, row.end_bp1), axis=1)
        df["in_viral_region2"] = np.where(
            df.contig_id1 == df.contig_id2,
            df.apply(lambda row: in_region(row.pos2_start, row.pos2_end, row.start_bp1, row.end_bp1), axis=1),
            0
        )

        df.to_csv(output.linkage_positions_annotated, sep='\t', index=False)

rule compile_linkages_between_viral_and_mags:
    input:
        linkage_positions_annotated = "results/{sample}/{sample}_linkage_positions_of_viral_mag_annot.tsv"
    output:
        virus_mag_linkages = "results/{sample}/{sample}_virus-mag_linkages.tsv"
    run:
        df = pd.read_csv(input.linkage_positions_annotated, sep='\t')

        # exclude viral self-linkages but include viral-nonviral linkages on the same contig
        df1 = df.loc[df['in_viral_region1'] + df['in_viral_region2'] == 1] \
            .groupby(['vir_id1', 'mag_id2'], as_index=False).agg(
                num_linkages = ('in_viral_region1', 'count')
            ).sort_values('num_linkages', ascending=False)

        df1.to_csv(output.virus_mag_linkages, sep='\t', index=False)

rule compile_linkages_between_viral:
    input:
        linkage_positions_annotated = "results/{sample}/{sample}_linkage_positions_of_viral_mag_annot.tsv"
    output:
        viral_HiC_stats = "results/{sample}/{sample}_viral_HiC_stats.tsv"
    run:
        df = pd.read_csv(input.linkage_positions_annotated, sep='\t')

        df['self_link'] = (df['in_viral_region1'] + df['in_viral_region2'] == 2).astype(int)
        df['nonself_link'] = (df['in_viral_region1'] + df['in_viral_region2'] == 1).astype(int)
        df['any_link'] = (df['in_viral_region1'] + df['in_viral_region2'] >= 1).astype(int)

        df1 = df.groupby('vir_id1', as_index=False).agg(
            Nself = ('self_link', 'sum'), M = ('nonself_link', 'sum'), N = ('any_link', 'sum')
        ).rename(columns={'vir_id1': 'vir_id'})

        df1.to_csv(output.viral_HiC_stats, sep='\t', index=False)

rule compute_virus_mag_linkages_pval:
    input:
        virus_mag_linkages = "results/{sample}/{sample}_virus-mag_linkages.tsv",
        viral_HiC_stats = "results/{sample}/{sample}_viral_HiC_stats.tsv",
        mag_HiC_stats = "results/{sample}/{sample}_mag_HiC_stats.tsv"
    output:
        virus_mag_linkages_pval = "results/{sample}/{sample}_virus-mag_linkages_pval.tsv"
    run:
        df = pd.read_csv(input.virus_mag_linkages, sep='\t')
        df_hic_v = pd.read_csv(input.viral_HiC_stats, sep='\t')
        df_hic_m = pd.read_csv(input.mag_HiC_stats, sep='\t')

        df = df.merge(df_hic_v.add_suffix('1'), on='vir_id1').merge(df_hic_m.add_suffix('2'), on='mag_id2')

        df['pval'] = df.apply(lambda row: binomial_pval(row.M1, row.p2, row.num_linkages), axis=1)

        pvals = df['pval'].values

        reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(pvals, method='fdr_bh')

        df['pval_corrected'] = pvals_corrected

        df.to_csv(output.virus_mag_linkages_pval, sep='\t', index=False)

# --- Aggregate results ---

rule combine_results:
    input:
        expand("results/{sample}/{sample}_virus-mag_linkages_pval.tsv", sample=SAMPLES)
    output:
        "results/virus-mag_HiC.tsv"
    run:
        dfs = []
        for s, f in zip(SAMPLES, input):
            df = pd.read_csv(f, sep='\t')
            df.insert(0, 'sample_id', s)
            dfs.append(df)

        pd.concat(dfs).to_csv(output[0], sep='\t', index=False)
