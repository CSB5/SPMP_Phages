"""
================================================================================
File:        bam_to_numlinkages.py
Description: Extracts number of intra and inter-contig linkages from a BAM file.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/HiC-host-association
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python bam_to_numlinkages.py -i <HiC.sorted.bam> -o <numlinkages.tsv>

Output:
    - TSV with columns: contig_id1, contig_id2, num_linkages

Notes:
    - This script *requires* a name-sorted BAM as input ('samtools sort -n').
    - Contig IDs within each pair are sorted (contig_id1 < contig_id2).
    - Contig pairs are sorted by decreasing number of linkages.
================================================================================
"""

import argparse
import sys
import pysam
from collections import Counter

def get_linkage_counts(bam_file: str) -> Counter:
    counts = Counter()
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        last_read = None
        for read in bam:
            if (not read.is_paired or
                read.is_unmapped or
                read.mate_is_unmapped or
                read.is_secondary or
                read.is_supplementary or
                read.is_duplicate):
                continue

            if last_read is None or read.query_name != last_read.query_name:
                last_read = read
                continue

            contig1 = last_read.reference_name
            contig2 = read.reference_name

            pair = tuple(sorted([contig1, contig2]))
            counts[pair] += 1
            last_read = None
    return counts

def main():
    parser = argparse.ArgumentParser(description="Extract number of linkages between contigs from a BAM file.")
    parser.add_argument("-i", "--input", default="-", help="Input BAM file (default: stdin)")
    parser.add_argument("-o", "--output", default="-", help="Output TSV file (default: stdout)")
    args = parser.parse_args()

    linkage_counts = get_linkage_counts(args.input)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write("contig_id1\tcontig_id2\tnum_linkages\n")
    for (c1, c2), n in linkage_counts.most_common():
        out.write(f"{c1}\t{c2}\t{n}\n")
    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()
