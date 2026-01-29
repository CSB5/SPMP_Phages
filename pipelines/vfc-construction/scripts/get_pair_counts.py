"""
================================================================================
File:        get_pair_counts.py
Description: Calculates pairwise similarities between viral genomes given
             protein cluster assignments.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/vfc-construction
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python get_pair_counts.py --mcl <pcs.mcl> --prot <proteins.list> --n_pc <num_pc> --out_pc <vr95_c50-PCs.tsv> --out_sim <vr95_c50-pairs.tsv>

Notes:
    - Pairwise similarities are defined as the mean % shared distinct protein
      clusters (PCs) between each pair of genomes.
    - Singleton PCs and orphans (proteins not in DIAMOND blastp results) are
      included in the total PC counts for each genome.
================================================================================
"""

import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Calculate pairwise similarities between viral genomes given protein clusters.')
    parser.add_argument('--mcl', required=True, help='MCL protein clustering output (.mcl)')
    parser.add_argument('--prot', required=True, help='List of all protein names (.list)')
    parser.add_argument('--n_pc', type=int, default=1, help='Minimum shared PCs to define a linkage between genomes (default: 1)')
    parser.add_argument('--out_pc', required=True, help='Table of PC assignments of all proteins.')
    parser.add_argument('--out_sim', required=True, help='Table of pairwise similarities between genomes')
    return parser.parse_args()

def get_df_pc(mclpcfile, protnamefile):
    print(f"[*] Reading protein clusters from {mclpcfile}...")
    
    all_rows = []
    with open(mclpcfile, 'r') as f:
        for i, line in enumerate(f):
            pc_id = f'PC_{i:06}'
            proteins = line.strip().split('\t')
            for prot in proteins:
                all_rows.append({'protein_id': prot, 'cluster': pc_id})
    
    df_clustered = pd.DataFrame(all_rows)
    
    print(f"[*] Reading {protnamefile} to retrieve orphan proteins...")
    df_all = pd.read_csv(protnamefile, names=['protein_id'])
    
    df_pc = df_all.merge(df_clustered, on='protein_id', how='left')
    
    mask = df_pc['cluster'].isnull()
    orphan_count = mask.sum()
    if orphan_count > 0:
        start_idx = len(df_clustered)
        orphan_ids = [f'PC_{i:06}' for i in range(start_idx, start_idx + orphan_count)]
        df_pc.loc[mask, 'cluster'] = orphan_ids

    df_pc['contig_id'] = df_pc['protein_id'].str.rsplit('_', n=1).str[0]
    
    return df_pc

def get_perc_shared_pcs(df_pc, n_pc_threshold):
    print("[*] Calculating % shared protein clusters...")
    
    # 1. Filter for clusters that appear in more than one contig
    cluster_counts = df_pc.groupby('cluster')['contig_id'].nunique()
    valid_clusters = cluster_counts[cluster_counts > 1].index
    
    df_filtered = df_pc[df_pc['cluster'].isin(valid_clusters)][['contig_id', 'cluster']].drop_duplicates()

    # 2. Perform a self-merge to find pairs sharing clusters
    df_pairs = pd.merge(df_filtered, df_filtered, on='cluster')
    
    # 3. Keep only unique pairs (A < B) to avoid self-links and duplicates (A-B and B-A)
    df_pairs = df_pairs[df_pairs['contig_id_x'] < df_pairs['contig_id_y']]
    
    # 4. Count shared clusters
    df_counts = df_pairs.groupby(['contig_id_x', 'contig_id_y']).size().reset_index(name='num_shared_pcs')
    
    # 5. Apply threshold
    df_counts = df_counts[df_counts['num_shared_pcs'] >= n_pc_threshold]
    
    # 6. Total distinct PCs per contig
    contig_pc_counts = df_pc.groupby('contig_id')['cluster'].nunique()
    
    df_counts = df_counts.merge(contig_pc_counts.rename('num_pc1'), left_on='contig_id_x', right_index=True)
    df_counts = df_counts.merge(contig_pc_counts.rename('num_pc2'), left_on='contig_id_y', right_index=True)
    
    # 7. Calculate percentages
    df_counts['perc_shared_pc1'] = (df_counts['num_shared_pcs'] / df_counts['num_pc1']) * 100
    df_counts['perc_shared_pc2'] = (df_counts['num_shared_pcs'] / df_counts['num_pc2']) * 100
    df_counts['mean_perc_shared_pcs'] = (df_counts['perc_shared_pc1'] + df_counts['perc_shared_pc2']) / 2
    
    return df_counts

def main():
    args = get_args()
    
    df_pc = get_df_pc(args.mcl, args.prot)
    df_pc.to_csv(args.out_pc, sep='\t', index=False)
    print(f"[+] Saved all protein clusters (including singletons and orphans) to {args.out_pc}")
    
    df_results = get_perc_shared_pcs(df_pc, args.n_pc)    
    df_results.rename(columns={'contig_id_x': 'vir_id1', 'contig_id_y': 'vir_id2'}).to_csv(args.out_sim, sep='\t', index=False)
    print(f"[+] Results saved to {args.out_sim}")

if __name__ == "__main__":
    main()
