"""
================================================================================
File:        vfc.smk
Description: Construction of viral clusters from a set of input and reference
             vOTUs based on protein sharing.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/vfc-construction
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    snakemake -s vfc.smk --cores 48

Inputs:
    - data/v95.fna
    - data/ref_fams95.fna

Notes:
    - Proteins are clustered using MCL with ANI*AF as similarity measure.
    - >=90% complete vOTUs are hierarchically clustered using UPGMA with
      (100 - mean % shared distinct PCs) as distance measure.
    - 50-90% complete vOTUs are recruited to VCs at d=80 (20% shared PCs).
================================================================================
"""

import pandas as pd

GENOMAD_DB = "resources/genomad_db"

THREADS_L = 48
THREADS_M = 24

DATASETS = ["v95", "ref_fams95"]

rule all:
    input:
        "results/vr95_c50-upgma.tsv"

# --- Completeness estimation ---

rule checkv:
    input: "data/{dataset}.fna"
    output:
        dir = directory("results/checkv_{dataset}"),
        table = "results/checkv_{dataset}/quality_summary.tsv"
    threads: THREADS_M
    shell:
        "checkv end_to_end {input} {output.dir} -t {threads}"

# --- Gene identification of >=50% complete vOTUs ---

rule filter_completeness_50:
    input: "results/checkv_{dataset}/quality_summary.tsv"
    output: "results/{dataset}_c50.list"
    run:
        df = pd.read_csv(input[0], sep='\t')
        df.loc[df.completeness >= 50][['contig_id']].to_csv(output[0], header=False, index=False)

rule get_fasta_from_list:
    input:
        list = "results/{dataset}_c50.list",
        fasta = "data/{dataset}.fna"
    output: "results/{dataset}_c50.fna"
    shell:
        "seqtk subseq {input.fasta} {input.list} > {output}"

rule genomad:
    input: "results/{dataset}_c50.fna"
    output:
        dir = directory("results/{dataset}_c50_ganno"),
        faa = "results/{dataset}_c50_ganno/{dataset}_c50_annotate/{dataset}_c50_proteins.faa"
    threads: THREADS_M
    shell:
        "genomad annotate -t {threads} {input} {output.dir} {GENOMAD_DB}"

# --- Generation of protein clusters ---

rule combine_proteins:
    input:
        proteins = expand("results/{ds}_c50_ganno/{ds}_c50_annotate/{ds}_c50_proteins.faa", ds=DATASETS)
    output: "results/vr95_c50.faa"
    shell:
        "cat {input.proteins} > {output}"

rule diamond_blastp:
    input: "results/vr95_c50.faa"
    output: "results/vr95_c50-blastp.tsv"
    params:
        db = "results/vr95_c50"
    threads: THREADS_L
    shell:
        """
        diamond makedb --in {input} -d {params.db}
        diamond blastp -d {params.db} -q {input} -o {output} --more-sensitive -e 1e-5 \
        --max-target-seqs 10000 -f 6 qseqid sseqid pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qlen slen --header --threads {threads}
        """

rule mcl_clustering:
    input: "results/vr95_c50-blastp.tsv"
    output: "results/vr95_c50-blastp.mcl"
    threads: THREADS_L
    shell:
        """
        awk '{{OFS = "\\t"}} NR>3 {{print $1, $2, $3*$4}}' {input} > results/vr95_c50-blastp.abc
        mcl results/vr95_c50-blastp.abc --abc -I 2 -te {threads} -o {output}
        """

# --- Computation of pairwise similarities ---

rule get_all_proteins:
    input:
        faa = "results/vr95_c50.faa"
    output:
        list = "results/vr95_c50-prot.list"
    shell:
        "grep '>' {input.faa} | awk '{{print $1}}' | tr -d '>' > {output.list}"

rule get_shared_proteins:
    input:
        mcl = "results/vr95_c50-blastp.mcl",
        list = "results/vr95_c50-prot.list"
    output:
        pairs = "results/vr95_c50-pairs.tsv",
        pcs = "results/vr95_c50-PCs.tsv"
    shell:
        "python scripts/get_pair_counts.py --mcl {input.mcl} --prot {input.list} --n_pc 1 --out_pc {output.pcs} --out_sim {output.pairs}"

# --- Hierarchical clustering of >=90% complete vOTUs ---

rule filter_completeness_90:
    input: 
        tables = expand("results/checkv_{ds}/quality_summary.tsv", ds=DATASETS)
    output: "results/vr95_c90.list"
    run:
        c90_ids = []
        for table in input.tables:
            df = pd.read_csv(table, sep='\t')
            
            ids = df.loc[df.completeness >= 90, 'contig_id'].tolist()
            c90_ids.extend(ids)
        
        with open(output[0], 'w') as f:
            for item in c90_ids:
                f.write(item + "\n")

rule get_completeness_90_pairs:
    input:
        pairs = "results/vr95_c50-pairs.tsv",
        c90_list = "results/vr95_c90.list"
    output:
        pairs = "results/vr95_c90-pairs.tsv",
        abc = "results/vr95_c90-pairs.abc"
    run:
        df = pd.read_csv(input.pairs, sep='\t')
        c90 = pd.read_csv(input.c90_list, names=['v']).v.tolist()
        df1 = df.loc[df.vir_id1.isin(c90) & df.vir_id2.isin(c90)]
        df1.to_csv(output.pairs, sep='\t', index=False)

        df1[['vir_id1','vir_id2','mean_perc_shared_pcs']].to_csv(output.abc, sep='\t', index=False, header=False)

rule run_hierarchical_clustering:
    input:
        abc = "results/vr95_c90-pairs.abc",
        c90_list = "results/vr95_c90.list"
    output: "results/vr95_c90-upgma.tsv"
    shell: "python scripts/perform_hierarchical_clustering.py -i {input.abc} -l {input.c90_list} -o {output}"

# --- Recruitment of 50-90% complete vOTUs into VCs at d=80 ---

rule recruit_50_90:
    input:
        pairs_50 = "results/vr95_c50-pairs.tsv",
        c90_list = "results/vr95_c90.list",
        upgma_90 = "results/vr95_c90-upgma.tsv"
    output:
        "results/vr95_c50-upgma.tsv"
    run:
        df = pd.read_csv(input.pairs_50, sep='\t')
        c90 = pd.read_csv(input.c90_list, names=['v']).v.tolist()
        
        df1 = df.loc[df.vir_id1.isin(c90) & ~df.vir_id2.isin(c90)][['vir_id1','vir_id2','perc_shared_pc1']].rename(columns={'vir_id2': 'vir_id'})
        df2 = df.loc[~df.vir_id1.isin(c90) & df.vir_id2.isin(c90)][['vir_id2','vir_id1','perc_shared_pc2']].rename(columns={'vir_id2': 'vir_id1', 'vir_id1': 'vir_id', 'perc_shared_pc2': 'perc_shared_pc1'})
        df1 = pd.concat([df1,df2], ignore_index=True)
        
        df_c90vc = pd.read_csv(input.upgma_90, sep='\t')
        df1 = df1.merge(df_c90vc.rename(columns={'vir_id': 'vir_id1'}), on='vir_id1')
        
        df2 = df1.groupby(['vir_id','VC80'], as_index=False).perc_shared_pc1.mean()
        df3 = df2.loc[df2.perc_shared_pc1 >= 30].copy()

        df4 = df3.sort_values(['vir_id','perc_shared_pc1'], ascending=[True,False]).drop_duplicates('vir_id')
        drange1 = [80,85,90,95]
        df_out = df4[['vir_id','VC80']].merge(df1[['vir_id'] + ['VC' + str(d) for d in drange1]], on=['vir_id','VC80'], how='left').drop_duplicates(['vir_id'] + ['VC' + str(d) for d in drange1])
        
        df90 = pd.read_csv(input.upgma_90, sep='\t')
        pd.concat([df90, df_out]).to_csv(output[0], sep='\t', index=False)
