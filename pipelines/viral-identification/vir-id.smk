"""
================================================================================
File:        vir-id.smk
Description: Viral identification leveraging VirSorter2, CheckV, and geNomad
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/viral-identification
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    snakemake -s vir-id.smk --cores 48

Inputs:
    - Metagenomic assembly files (modify config.yaml and/or `SAMPLES` variable
      to match the sample names)
================================================================================
"""

configfile: "config/config.yaml"

# sample names
SAMPLES = [config["sample_prefix"] + f"{n:0{config['padding']}d}" for n in range(1, config["num_samples"] + 1)]

import pandas as pd
from Bio import SeqIO

rule all:
    input:
        "results/viral.fna",
        "results/v95.fna"

# --- VirSorter2 ---

rule vs2:
    input:
        config["contigs_dir"] + "{sample}" + config["contigs_filename"]
    output:
        dir = directory("results/vs2-checkv/{sample}/vs2"),
        fna = "results/vs2-checkv/{sample}/vs2/final-viral-combined.fa"
    params:
        l=config.get("min_len", 5000),
        s=config.get("min_vs2_score", 0.5)
    threads: 16
    shell:
        "virsorter run --keep-original-seq -i {input} -w {output.dir} "
        "--include-groups dsDNAphage,ssDNA --min-length {params.l} --min-score {params.s} -j {threads} all"

# --- CheckV ---

rule checkv:
    input:
        "results/vs2-checkv/{sample}/vs2/final-viral-combined.fa"
    output:
        dir = directory("results/vs2-checkv/{sample}/checkv"),
        pv = "results/vs2-checkv/{sample}/checkv/proviruses.fna",
        v = "results/vs2-checkv/{sample}/checkv/viruses.fna"
    threads: 12
    shell:
        "checkv end_to_end {input} {output.dir} -t {threads}"

rule get_checkv_out:
    input:
        "results/vs2-checkv/{sample}/checkv/proviruses.fna",
        "results/vs2-checkv/{sample}/checkv/viruses.fna"
    output:
        "results/vs2-checkv/{sample}/checkv_out.fna"
    shell:
        # merge and clean headers (removing text after space)
        "cat {input} | sed 's/ .*//' > {output}"

rule checkv_2:
    input:
        "results/vs2-checkv/{sample}/checkv_out.fna"
    output:
        dir = directory("results/vs2-checkv/{sample}/checkv_2"),
        tsv = "results/vs2-checkv/{sample}/checkv_2/quality_summary.tsv"
    threads: 4
    shell:
        "checkv end_to_end {input} {output.dir} -t {threads}"

# --- Filter length >= min_len and viral >= host genes ---

rule filter_and_rename:
    input:
        fna = "results/vs2-checkv/{sample}/checkv_out.fna",
        tsv = "results/vs2-checkv/{sample}/checkv_2/quality_summary.tsv"
    output:
        fna = "results/vs2-checkv/{sample}/filtered.fna",
        tsv = "results/vs2-checkv/{sample}/filtered-info.tsv"
    params:
        min_len = config.get("min_len", 5000),
        prefix_to_remove = config.get("contig_name_prefix", "")
    run:
        df = pd.read_csv(input.tsv, sep='\t')
        filtered_df = df[(df['contig_length'] >= params.min_len) & (df['viral_genes'] >= df['host_genes'])].copy()

        filtered_ids = set(filtered_df['contig_id'])
        records_to_save = []
        
        for record in SeqIO.parse(input.fna, "fasta"):
            if record.id in filtered_ids:
                new_id = f"{wildcards.sample}_{record.id}".replace(params.prefix_to_remove, "").replace("|", "_")
                record.id = new_id
                record.description = ""
                records_to_save.append(record)
        
        SeqIO.write(records_to_save, output.fna, "fasta")
        
        filtered_df['contig_id'] = filtered_df['contig_id'].apply(
            lambda x: f"{wildcards.sample}_{x}".replace(params.prefix_to_remove, "").replace("|", "_")
        )
        filtered_df.to_csv(output.tsv, sep='\t', index=False)

# --- geNomad ---

rule genomad:
    input:
        config["contigs_dir"] + "{sample}" + config["contigs_filename"]
    output:
        dir = directory("results/genomad/{sample}"),
        fna = "results/genomad/{sample}/{sample}_summary/{sample}_virus.fna",
        tsv = "results/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
    params:
        db = config.get("genomad_db", "")
    threads: 16
    shell:
        "genomad end-to-end -t {threads} --enable-score-calibration {input} {output.dir} {params.db}"

# --- Filter length >= min_len and FDR < max_gmd_fdr ---

rule filter_genomad:
    input:
        fna = "results/genomad/{sample}/{sample}_summary/{sample}_virus.fna",
        tsv = "results/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
    output:
        fna = "results/genomad/{sample}/gvirs.flt.fna",
        tsv = "results/genomad/{sample}/gvirs.flt.tsv"
    params:
        min_len = config.get("min_len", 5000),
        max_fdr = config.get("max_gmd_fdr", 0.01)
    run:
        df = pd.read_csv(input.tsv, sep='\t')
        filtered_df = df[(df['length'] >= params.min_len) & (df['fdr'] < params.max_fdr)].copy()

        filtered_ids = set(filtered_df['seq_name'])
        records_to_save = []
        
        for record in SeqIO.parse(input.fna, "fasta"):
            clean_id = record.id.split(' ')[0]
            
            if clean_id in filtered_ids:
                new_id = f"{wildcards.sample}_{clean_id}".replace("|", "_")
                record.id = new_id
                record.description = ""
                records_to_save.append(record)

        SeqIO.write(records_to_save, output.fna, "fasta")

        filtered_df['seq_name'] = filtered_df['seq_name'].apply(
            lambda x: f"{wildcards.sample}_{x}".replace("|", "_")
        )
        filtered_df.to_csv(output.tsv, sep='\t', index=False)

# --- Aggregate and dereplicate viral sequences ---

rule get_viral_sequences:
    input:
        genomad = expand("results/genomad/{sample}/gvirs.flt.fna", sample=SAMPLES),
        vs2_checkv = expand("results/vs2-checkv/{sample}/filtered.fna", sample=SAMPLES)
    output:
        "results/v.fna"
    shell:
        "cat {input} > {output}"

rule cdhit:
    input:
        "results/v.fna"
    output:
        viral = "results/viral.fna",
        v99 = "results/v99.fna",
        v95 = "results/v95.fna"
    threads: 24
    shell:
        """
        cd-hit-est -i {input} -o {output.viral} -c 1 -n 10 -d 0 -M 0 -T {threads}
        cd-hit-est -i {output.viral} -o {output.v99} -c 0.99 -n 10 -d 0 -M 0 -T {threads}
        cd-hit-est -i {output.v99} -o {output.v95} -c 0.95 -G 0 -aS 0.85 -n 10 -d 0 -M 0 -T {threads}
        """
