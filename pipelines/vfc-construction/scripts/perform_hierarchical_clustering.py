"""
================================================================================
File:        perform_hierarchical_clustering.py
Description: Performs UPGMA hierarchical clustering on a pairwise similarity
             matrix to generate cluster labels (VCs) at a range of distance
             thresholds.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/vfc-construction
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python perform_hierarchical_clustering.py -i <input.abc> -l <ids.list> -o <upgma.tsv>

Notes:
    - VC labels are not comparable across different thresholds, e.g., VC 50 at
      d=60 is not related to VC 50 at d=80.
================================================================================
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import average, fcluster
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Perform hierarchical clustering based on shared protein clusters.')
    parser.add_argument('-i', '--input', required=True, help='Input pairwise similarity file (.abc)')
    parser.add_argument('-l', '--list', required=True, help='List of all viral IDs (.list)')
    parser.add_argument('-o', '--output', required=True, help='Table of VC labels at multiple thresholds')
    parser.add_argument('-t', '--thresholds', nargs='+', type=float,
                        default=[60, 65, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 85, 90, 95],
                        help='Distance thresholds to cut the dendrogram')
    return parser.parse_args()

def main():
    args = get_args()

    print(f"[*] Loading viral ID list: {args.list}")
    df_v = pd.read_csv(args.list, names=['vir_id'])
    df_v['index'] = range(len(df_v))
    m = len(df_v)
    
    print(f"[*] Loading pairwise similarities: {args.input}")
    df_sim = pd.read_csv(args.input, sep='\t', names=['vir_id1', 'vir_id2', 'sim'])

    # 1. Prepare indices for the condensed matrix
    # Map IDs to their integer positions
    id_map = dict(zip(df_v['vir_id'], df_v['index']))
    df_sim['i'] = df_sim['vir_id1'].map(id_map)
    df_sim['j'] = df_sim['vir_id2'].map(id_map)

    # Ensure i < j for condensed matrix indexing
    mask = df_sim['i'] > df_sim['j']
    df_sim.loc[mask, ['i', 'j']] = df_sim.loc[mask, ['j', 'i']].values

    # 2. Construct Condensed Distance Matrix
    # Formula for index in condensed matrix: m*i + j - ((i+2)*(i+1))//2
    n = m * (m - 1) // 2
    print(f"[*] Initializing distance matrix for {m} viruses ({n} entries)...")
    
    D = np.full(n, 100.0, dtype=np.float32)

    i = df_sim['i'].values
    j = df_sim['j'].values
    entries = m * i + j - (i + 2) * (i + 1) // 2
    
    D[entries] = 100.0 - df_sim['sim'].values
    print("[+] Distance matrix constructed.")

    # 3. Perform UPGMA
    print("[*] Performing UPGMA clustering...")
    Z = average(D)
    
    # 4. Cut dendrogram at various thresholds
    results_df = df_v[['vir_id']].copy()
    
    print(f"[*] Cutting dendrogram at d={args.thresholds}...")
    for d in args.thresholds:
        # fcluster returns cluster IDs (1-based)
        cluster_labels = fcluster(Z, d, criterion='distance')
        results_df[f'VC{int(d)}'] = cluster_labels - 1

    # 5. Save results
    results_df.to_csv(args.output, sep='\t', index=False)
    print(f"[+] Clustering results saved to {args.output}")

if __name__ == "__main__":
    main()
