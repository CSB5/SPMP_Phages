"""
================================================================================
File:        solve_hic_model.py
Description: Solves Hi-C noise model for MAGs.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/HiC-host-association
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python solve_hic_model.py -i <maglinkages.tsv> -o <mag_HiC_stats.tsv>

Output:
    - TSV with columns: mag_id, Nself, M, max_othermag_m, N, e, p

================================================================================
"""

import pandas as pd
import numpy as np
import argparse

def solve_model(input_file, output_file, e_init, max_iter, tol):
    df = pd.read_csv(input_file, sep='\t')

    # get self-linkages
    df_sum = df.loc[df.mag_id1 == df.mag_id2].copy()
    df_sum = df_sum.rename(columns={'mag_id1': 'mag_id', 'num_linkages': 'Nself'})[['mag_id', 'Nself']].sort_values('Nself', ascending=False)

    # get nonself-linkages
    df1 = df.loc[df.mag_id1 != df.mag_id2]
    df_long = pd.melt(
        df1,
        id_vars='num_linkages',
        value_vars=['mag_id1', 'mag_id2'],
        value_name='mag_id'
    )
    df_sum1 = df_long.groupby('mag_id', as_index=False)['num_linkages'].agg(
        M='sum', 
        max_othermag_m='max'
    )
    
    df_sum = df_sum.merge(df_sum1, on='mag_id', how='left').fillna(0)

    # total linkages
    df_sum['N'] = df_sum.Nself + df_sum.M

    M = df_sum.M.tolist()
    N = df_sum.N.tolist()

    # optimization loop
    e = e_init * np.ones_like(M, dtype=float)

    converged = False
    for it in range(max_iter):
        e_old = e.copy()

        # Step 1: update p_s
        p = N * e / np.sum(N * e)

        # Step 2: update e_s
        e = M / (N * (1 - p))

        # check convergence
        if np.all(np.abs(e - e_old) < tol):
            print(f"Converged after {it+1} iterations.")
            converged = True
            break

    if not converged:
        print(f"Warning: Did not converge within {max_iter} iterations.")

    df_sum['e'] = e
    df_sum['p'] = p

    df_sum.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Solve Hi-C null model for MAGs.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV of inter-MAG linkages")
    parser.add_argument("-o", "--output", required=True, help="Output TSV of inferred model parameters")
    parser.add_argument("--e_init", type=float, default=0.01, help="Initial value of e")
    parser.add_argument("--max_iter", type=int, default=1000, help="Max no. of iterations")
    parser.add_argument("--tol", type=float, default=1e-8, help="Convergence tolerance")
    args = parser.parse_args()

    solve_model(args.input, args.output, args.e_init, args.max_iter, args.tol)

if __name__ == "__main__":
    main()
