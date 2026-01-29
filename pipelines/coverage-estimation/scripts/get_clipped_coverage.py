"""
================================================================================
File:        get_clipped_coverage.py
Description: Calculates clipped mean and standard deviation of base coverage
             for a set of references in a sample.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/coverage-estimation
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python get_clipped_coverage.py -c <coverage.tsv> -d <depth.tsv> -s <sample_id> -o <clippedcov.tsv> -p <prefix>

Notes:
    - Zero-padding of base coverage values is performed if samtools depth was
      not run with -a option.
================================================================================
"""

import argparse
import pandas as pd
import numpy as np
from typing import List

def calculate_clipped_stats(coverage_values: List[int], total_len: int, threshold: float, depth_all: bool):
    padded_values = coverage_values.copy()
    if not depth_all:
        padded_values += [0] * (total_len - len(coverage_values))

    padded_values.sort()

    frac = threshold / 100
    clip_l = int(np.ceil(frac * total_len))
    clip_h = int(np.floor((1 - frac) * total_len))

    clipped_seq = padded_values[clip_l:clip_h]

    if not clipped_seq:
        return 0.0, 0.0

    return np.mean(clipped_seq), np.std(clipped_seq)

def main():
    parser = argparse.ArgumentParser(description='Calculate clipped coverage statistics.')
    parser.add_argument('-c', '--coverage_file', required=True, help='samtools coverage file')
    parser.add_argument('-d', '--depth_file', required=True, help='samtools depth file')
    parser.add_argument('-s', '--sample_id', required=True, help='Sample ID')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-t', '--clip_threshold', type=float, default=10, help='Top and bottom percentile base coverage values to clip (default: 10)')
    parser.add_argument('-a', '--depth_all', action="store_true", help='If samtools depth was run with -a option')
    parser.add_argument('-p', '--header_prefix', default='contig', help='Header prefix for reference id and len (default: contig)')
    args = parser.parse_args()

    df_cov = pd.read_csv(args.coverage_file, sep='\t')
    df_dep = pd.read_csv(args.depth_file, sep='\t', names=['#rname', 'pos', 'depth'])

    results = []

    for ref, df in df_dep.groupby('#rname', sort=False):
        cov_row = df_cov[df_cov['#rname'] == ref].iloc[0]
        ref_len = int(cov_row['endpos'])
        num_reads = int(cov_row['numreads'])
        covered_bases = int(cov_row['covbases'])
        cov_breadth = float(cov_row['coverage'])

        if not args.depth_all and len(df) != covered_bases:
            print(f"Warning: samtools coverage and depth have mismatched coverage breadth for {ref} in sample {args.sample_id}.")
        elif args.depth_all and len(df) != ref_len:
            print(f"Warning: samtools coverage and depth have mismatched reference length for {ref} in sample {args.sample_id}.")

        depth_list = df['depth'].to_list()
        c_mean, c_std = calculate_clipped_stats(depth_list, ref_len, args.clip_threshold, args.depth_all)

        results.append({
            f'{args.header_prefix}_id': ref,
            f'{args.header_prefix}_len': ref_len,
            'sample_id': args.sample_id,
            'numreads': num_reads,
            'covbases': covered_bases,
            'covbreadth': cov_breadth,
            'covmean': c_mean,
            'covstd': c_std
        })

    output_df = pd.DataFrame(results)
    output_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
