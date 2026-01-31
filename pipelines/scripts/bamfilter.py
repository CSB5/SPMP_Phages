"""
================================================================================
File:        bamfilter.py
Description: Filters BAM reads by minimum alignment identity, coverage, and
             proper pairing.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/coverage-estimation
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python bamfilter.py -i <in.bam> -o <out.bam> --min_id <id> --min_cov <cov> --check_proper_pair <1>

================================================================================
"""

import argparse
import pysam

def is_valid_read(read: pysam.AlignedSegment, min_cov: float, min_id: float, paired_flag: bool) -> bool:
    if not read.is_mapped:
        return False
    if paired_flag and not read.is_proper_pair:
        return False

    len_ops, num_ops = read.get_cigar_stats()
    read_length = read.infer_read_length()
    if read_length == 0:
        return False

    coverage = sum(len_ops[0:2]) / read_length * 100        # ([alignment length] + [insertion length]) / [read length] * 100
    identity = (1 - len_ops[10] / sum(len_ops[0:3])) * 100  # (1 - [edit distance] / ([alignment length] + [insertion length] + [deletion length])) * 100

    return coverage >= min_cov and identity >= min_id

def filter_bam(input_bam: str, output_bam: str, min_cov: float, min_id: float, paired_flag: bool):
    with pysam.AlignmentFile(input_bam, "rb") as bam_in, \
         pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out:
        for read in bam_in:
            if is_valid_read(read, min_cov, min_id, paired_flag):
                bam_out.write(read)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM reads by minimum alignment identity, coverage, and proper pairing.")
    parser.add_argument("-i", default="-", help="Input BAM file (default: stdin)")
    parser.add_argument("-o", default="-", help="Output BAM file (default: stdout)")
    parser.add_argument("--min_id", type=float, default=95, help="Minimum alignment identity [0,100] (default: 95)")
    parser.add_argument("--min_cov", type=float, default=80, help="Minimum read coverage [0,100] (default: 80)")
    parser.add_argument("--check_proper_pair", type=int, default=1, help="Check if reads are properly paired (default: 1)")
    args = parser.parse_args()

    paired_flag = args.check_proper_pair == 1
    filter_bam(args.i, args.o, args.min_cov, args.min_id, paired_flag)

if __name__ == "__main__":
    main()
