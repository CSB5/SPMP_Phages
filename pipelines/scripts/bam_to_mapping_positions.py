"""
================================================================================
File:        bam_to_mapping_positions.py
Description: Extracts mapping positions between pairs of contigs belonging to
             two specified lists from a BAM file.
Author:      Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/HiC-host-association
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python bam_to_mapping_positions.py -i <HiC.sorted.bam> -l1 <contigs1.list> -l2 <contigs2.list> -o <mapping_positions.tsv>

Output:
    - TSV with columns: contig_id1, pos1_start, pos1_end, contig_id2, pos2_start, pos2_end

Notes:
    - This script *requires* a name-sorted BAM as input ('samtools sort -n').
================================================================================
"""

import argparse
import sys
import pysam

def get_mapping_positions(bam_file, contigs1, contigs2):
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

            r1, r2 = (last_read, read) if last_read.is_read1 else (read, last_read)

            contig1 = r1.reference_name
            contig2 = r2.reference_name

            pos1_start = r1.reference_start + 1
            pos1_end = r1.reference_end
            pos2_start = r2.reference_start + 1
            pos2_end = r2.reference_end

            if contig1 in contigs1 and contig2 in contigs2:
                yield [contig1, pos1_start, pos1_end, contig2, pos2_start, pos2_end]
            elif contig1 in contigs2 and contig2 in contigs1:
                yield [contig2, pos2_start, pos2_end, contig1, pos1_start, pos1_end]

            last_read = None

def main():
    parser = argparse.ArgumentParser(description="Extract mapping positions between contigs from a BAM file, optionally mapping to a specified list of contigs.")
    parser.add_argument("-i", "--input", default="-", help="Input BAM file (default: stdin)")
    parser.add_argument("-l1", "--list1", help="Reference contig list 1. Read pairs where one read maps to a contig in this list are extracted")
    parser.add_argument("-l2", "--list2", help="Reference contig list 2. Read pairs where one read maps to a contig in this list are extracted")
    parser.add_argument("-o", "--output", default="-", help="Output TSV file (default: stdout)")
    args = parser.parse_args()

    contigs1 = set(line.strip() for line in open(args.list1))
    contigs2 = set(line.strip() for line in open(args.list2))

    mapping_positions = get_mapping_positions(args.input, contigs1, contigs2)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write("contig_id1\tpos1_start\tpos1_end\tcontig_id2\tpos2_start\tpos2_end\n")
    for pos in get_mapping_positions(args.input, contigs1, contigs2):
        out.write("\t".join(map(str, pos)) + "\n")
    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()
