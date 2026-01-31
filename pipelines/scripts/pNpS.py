"""
================================================================================
File:        pNpS.py
Description: Process coverage and VCF files for a sample to produce SNP, codon,
             and pN/pS statistics.
Authors:     Windsor Koh, Hanrong Chen (chenhr@a-star.edu.sg)
Repository:  https://github.com/CSB5/SPMP_Phages/pipelines/pNpS
Citation:    See README.md for citation instructions.
================================================================================
Usage:
    python pNpS.py --sample <sample_id> --gene_file <processed_genes.tsv> \
                   --fasta_file <ref.fna> --cov_file <coverage.tsv> --cov_threshold <thresh> --vcf_file <in.vcf> \
                   --vcf_out <annotated.vcf> --gene_pNpS <gene_pNpS.tsv>

Prerequisites:
    - processed_genes.tsv should have been preprocessed using preprocess_gene_file_possible_sites().
    - in.vcf should have 1 row per ALT allele at each site (comma-separated ALT fields should be preprocessed).

Notes:
    - Codon statistics are computed for each gene a SNP is in.
    - Only sites with >1 allele are considered. Those with AF = 1 are discarded.
    - Synonymous/nonsynonymous mutations are checked w.r.t. the reference allele, or the first ALT allele if the reference AF = 0.
    - Genetic codes 11, 4, and 15 can be used.
    - Pseudocounts are added to the number of observed S and NS before calculating pN/pS.
    - Stop codons are excluded from the pN/pS calculation.
    - All sequences detected in the sample (coverage >= 70%) are returned, even if they have 0 SNPs.
    - codon_starting_index is w.r.t. the forward strand (even when the gene is on the reverse strand). Codons are obtained by reverse_complement(seq[i:i+3]).
================================================================================
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import os

# --- Process samtools coverage output to retain reference sequences satisfying coverage breadth threshold ---

def process_coverage_file(cov_file, covthresh=70):
    cov_stats = pd.read_csv(cov_file, sep = '\t').rename(columns={'#rname': 'CHROM'})
    cov_stats = cov_stats[cov_stats.coverage >= covthresh]
    cov_stats.set_index('CHROM', inplace=True)
    return cov_stats

# --- Generator of FASTA header and sequence ---

def read_fasta(fasta_file):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        yield record.id, str(record.seq.upper()) # convert to uppercase

# --- Process VCF file to retain SNPs belonging to references satisfying coverage breadth threshold, extract AF, keep AF < 1 ---

def process_vcf_file(vcf_file, cov_stats):
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

    vcf_df_filtered = vcf_df[vcf_df['CHROM'].isin(cov_stats.index)]

    # assumes VCF file has multiple rows for the same POS if there is >1 ALT allele (NOT comma-separated), or else returns only the first one
    vcf_df_filtered["AF"] = vcf_df_filtered["INFO"].str.extract(r"AF=([0-9\.eE+-]+)").astype(float)

    # keep AF < 1
    vcf_df_filtered = vcf_df_filtered[vcf_df_filtered["AF"].lt(1)]

    vcf_df_filtered['gene'] = ''
    vcf_df_filtered['codon_starting_index'] = 0
    vcf_df_filtered['snp_position_in_codon'] = 0
    vcf_df_filtered['ref_codon'] = ''
    vcf_df_filtered['mut_codon'] = ''
    vcf_df_filtered['mutation_type'] = ''
    return vcf_df_filtered

# --- Codon helper functions ---

gencode = {11:
{
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}, 15:
{
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'Q', # TAG: * -> Q
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}, 4:
{
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'W', 'TGG':'W' # TGA: * -> W
}}

# returns 0-based index for starting position of codon containing SNP
def get_codon_starting_index(snp_position, gene_start):
    codon_start_offset = (snp_position - gene_start) // 3 * 3
    return gene_start + codon_start_offset - 1

def mutate_codon(codon, pos, m):
    return codon[:pos] + m + codon[pos+1:]

def complement(m):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return comp.get(m, 'N') # returns 'N' as complement of non-standard characters

def reverse_complement(codon):
    return ''.join(complement(b) for b in codon[::-1])

def comp_codons(c1, c2, genetic_code=11):
        c1_aa = gencode[genetic_code].get(c1, None)
        c2_aa = gencode[genetic_code].get(c2, None)
        if c1_aa is None or c2_aa is None:
            return 'Unknown' # non-standard codon
        return 'Syn' if c1_aa == c2_aa else 'Non-syn'

# --- Return reference and mutated codons given by a string of concatenated SNPs ---

def get_ref_mutated_codons(snp_position, gene_start, strand, ref_sequence, SNPs=''):
    codons = []
    codon_starting_index = get_codon_starting_index(snp_position, gene_start)
    if strand == 1:
        ref_codon = ref_sequence[codon_starting_index:codon_starting_index+3]
        codons.append(ref_codon)
        snp_position_in_codon = (snp_position - gene_start)%3
        for m in SNPs:
            codons.append(mutate_codon(ref_codon, snp_position_in_codon, m))
    elif strand == -1:
        ref_codon = reverse_complement(ref_sequence[codon_starting_index:codon_starting_index+3])
        codons.append(ref_codon)
        snp_position_in_codon = 2 - (snp_position - gene_start)%3
        for m in SNPs:
            codons.append(mutate_codon(ref_codon, snp_position_in_codon, complement(m)))
    return codons, codon_starting_index, snp_position_in_codon

# --- Generate codon information of SNPs in a reference ---

def annotate_codons_of_ref(vcf_ref, ref_genes, ref_sequence):
    vcf_out = []

    for pos, vcf_pos in vcf_ref.groupby('POS'):
        relevant_genes = ref_genes[(ref_genes['start'] <= pos) & (ref_genes['end'] >= pos)]
        if len(relevant_genes) == 0: # SNP is intergenic
            vcf_pos.loc[vcf_pos['POS'] == pos, 'mutation_type'] = 'Intergenic'
            vcf_out.append(vcf_pos)
        else:
            for index, row in relevant_genes.iterrows(): # loop through every gene that SNP is found in
                vcf_gene = vcf_pos.copy()
                vcf_gene['gene'] = row['gene']
                SNPs = ''.join(vcf_gene['ALT'])
                codons, codon_starting_index, snp_position_in_codon = get_ref_mutated_codons(pos, row['start'], row['strand'], ref_sequence, SNPs)
                vcf_gene['codon_starting_index'] = codon_starting_index
                vcf_gene['snp_position_in_codon'] = snp_position_in_codon
                if (row['strand'] == 1 and codon_starting_index+3 == row['end']) or (row['strand'] == -1 and codon_starting_index+1 == row['start']): # SNP is in last codon of gene
                    vcf_gene['mutation_type'] = 'Stop'
                else:
                    if vcf_gene['AF'].sum() == 1: # AF of ALT alleles sum to 1, use first ALT allele as reference instead
                        vcf_gene = vcf_gene.iloc[1:]
                        codons = codons[1:]
                    vcf_gene['ref_codon'] = codons[0]
                    vcf_gene['mut_codon'] = codons[1:] # LHS and RHS should have the same length
                    vcf_gene['mutation_type'] = vcf_gene.apply(lambda r: comp_codons(r['ref_codon'], r['mut_codon'], row['genetic_code']), axis=1)
                vcf_out.append(vcf_gene)

    return pd.concat(vcf_out, ignore_index=True)

# --- Generate codon information of all SNPs in VCF ---

def annotate_codons(vcf, genes, fasta_file):
    vcf_refs = []
    for ref_name, ref_sequence in read_fasta(fasta_file):
        vcf_ref = vcf[vcf['CHROM'] == ref_name]
        if not vcf_ref.empty:
            genes_ref = genes[genes['CHROM'] == ref_name]
            vcf_ref = annotate_codons_of_ref(vcf_ref, genes_ref, ref_sequence)
            vcf_refs.append(vcf_ref)
    return pd.concat(vcf_refs, ignore_index=True)

# --- Count the number of synonymous and nonsynonymous single-base substitutions in a codon ---

def codon_possible_sites(codon, genetic_code=11):
    aa_orig = gencode[genetic_code].get(codon, None)
    if aa_orig is None: # if codon contains non-standard characters, skip it
        return 0, 0

    S = 0
    NS = 0
    for i in range(3):
        for m in 'ATCG':
            new_codon = mutate_codon(codon, i, m)
            if new_codon != codon: # only check 3x3 mutations
                if gencode[genetic_code][new_codon] == aa_orig:
                    S += 1
                else:
                    NS += 1
    return S, NS

# --- Find the number of possible S and NS mutations in a gene, skipping final stop codon ---

def gene_possible_sites(ref_sequence, gene_start, gene_end, strand, genetic_code=11):
    S_gene = 0
    NS_gene = 0
    if strand == 1:
        for i in range(gene_start-1, gene_end-3, 3): # loop through codons, skip last one
            S, NS = codon_possible_sites(ref_sequence[i:i+3], genetic_code)
            S_gene += S
            NS_gene += NS
    elif strand == -1:
        for i in range(gene_end-3, gene_start-1, -3): # loop through codons, skip last one
            S, NS = codon_possible_sites(reverse_complement(ref_sequence[i:i+3]), genetic_code)
            S_gene += S
            NS_gene += NS
    return S_gene, NS_gene

# --- Create 'CHROM', 'possible_S', 'possible_NS' columns in gene file ---

def preprocess_gene_file_possible_sites(gene_file, fasta_file):
    gene_table = pd.read_csv(gene_file, sep='\t')
    gene_table.insert(0, 'CHROM', gene_table['gene'].str.rsplit('_', n=1).str[0])

    fasta_sequences = {name: seq for name, seq in read_fasta(fasta_file)}

    def process_row(row):
        seq = fasta_sequences.get(row.CHROM, "")
        return gene_possible_sites(seq, row.start, row.end, row.strand, row.genetic_code)

    results = [process_row(row) for row in gene_table.itertuples(index=False)]
    gene_table[['possible_S', 'possible_NS']] = pd.DataFrame(results, index=gene_table.index)

    return gene_table

# --- Computed observed number of S and NS mutations, and pN/pS, of gene ---

def get_gene_pnps(gene_table, processed_vcf, cov_stats):
    gene_table_flt = gene_table[gene_table.CHROM.isin(cov_stats.index)]

    genes_all = []
    for ref_name, genes_ref in gene_table_flt.groupby('CHROM', sort=False):
        vcf_ref = processed_vcf[processed_vcf['CHROM'] == ref_name]
        for index, row in genes_ref.iterrows():
            vcf_gene = vcf_ref[vcf_ref['gene'] == row['gene']]
            genes_ref.at[index, 'num_SNP_pos'] = vcf_gene['POS'].nunique() # number of SNP positions
            genes_ref.at[index, 'S'] = len(vcf_gene[vcf_gene['mutation_type'] == 'Syn'])
            genes_ref.at[index, 'NS'] = len(vcf_gene[vcf_gene['mutation_type'] == 'Non-syn'])
        genes_all.append(genes_ref)

    genes_merged = pd.concat(genes_all, ignore_index=True)
    genes_merged['S_pseudocount'] = genes_merged['S'] + 1
    genes_merged['NS_pseudocount'] = genes_merged['NS'] + 1
    genes_merged['pNpS'] = (genes_merged['NS_pseudocount'] / genes_merged['possible_NS']) / (genes_merged['S_pseudocount'] / genes_merged['possible_S'])

    columns_to_int = ['num_SNP_pos', 'S', 'NS', 'S_pseudocount', 'NS_pseudocount']
    genes_merged[columns_to_int] = genes_merged[columns_to_int].astype(int)
    return genes_merged

def main():
    parser = argparse.ArgumentParser(description='Compute codon and pN/pS statistics.')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--cov_file', help='samtools coverage output')
    parser.add_argument('--vcf_file', help='VCF file')
    parser.add_argument('--gene_file', help='Preprocessed gene annotations')
    parser.add_argument('--fasta_file', help='Reference FASTA')
    parser.add_argument('--cov_threshold', type=float, default=70, help='Coverage threshold for detection in a sample (default: 70)')
    parser.add_argument('--vcf_out', help='Processed VCF file')
    parser.add_argument('--gene_pNpS', help='Gene-level pN/pS statistics')
    args = parser.parse_args()

    processed_genes = pd.read_csv(args.gene_file, sep='\t')
    cov_stats = process_coverage_file(args.cov_file, covthresh=args.cov_threshold)
    vcf = process_vcf_file(args.vcf_file, cov_stats)

    processed_vcf = annotate_codons(vcf, processed_genes, args.fasta_file)
    gene_pNpS = get_gene_pnps(processed_genes, processed_vcf, cov_stats)

    processed_vcf.insert(0, 'sample_id', args.sample)
    gene_pNpS.insert(0, 'sample_id', args.sample)

    processed_vcf.to_csv(args.vcf_out, sep='\t', index=False)
    gene_pNpS.to_csv(args.gene_pNpS, sep='\t', index=False)

if __name__ == "__main__":
    main()
