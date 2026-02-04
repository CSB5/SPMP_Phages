# Variant calling and pN/pS analysis

This is a bioinformatics pipeline for variant calling and pN/pS analysis of vOTUs. Genes of input vOTUs are identified, and the number of possible synonymous and nonsynonymous single-nucleotide substitutions is calculated for each gene. Paired-end Illumina reads are mapped to a database comprising input vOTUs and *decoy MAGs with prophages removed*. Only reads mapping to vOTUs are retained and filtered. Variants are called using LoFreq. Positions with >1 variants are retained, classified as intergenic, synonymous or nonsynonymous, and pN/pS of all genes in detected vOTUs are computed.

## 1. Overview

This pipeline performs the following steps:

1. **Gene annotation:** Genes are annotated using geNomad.
2. **Synonymous/nonsynonymous site counts:** The number of possible synonymous and nonsynonymous single-nucleotide substitutions is computed for each gene.
3. **Read mapping:** Paired-end Illumina reads are mapped to input vOTUs and decoy MAGs with prophages removed using BWA-MEM.
4. **BAM filtering:** Only reads mapping to input vOTUs are retained, and are filtered based on mapping identity (default: 95%), coverage (default: 80%), and proper pairing.
5. **Coverage computation:** The coverage breadth of each sequence is computed. vOTUs passing a coverage breadth threshold (default: 70%) are considered present.
6. **Variant calling:** Variants are called using LoFreq.
7. **Classification of variants:** Positions with >1 variants are classified as intergenic, synonymous or nonsynonymous.
8. **pN/pS computation:** pN/pS of all genes in detected vOTUs are computed.
9. **Aggregation of results:** Results are aggregated across samples.

## 2. Directory structure

```text
pNpS/
├── README.md
├── pNpS.smk
├── config.yaml
├── samples.tsv     # table of sample IDs and read FASTQ path names
└── results/        # outputs (not included)
data/
├── v95.fna         # input vOTUs (not included)
├── decoy_mags.fna  # decoy MAGs with prophages removed (not included)
└── reads/          # directory containing Illumina reads (not included)
scripts/
├── bamfilter.py    # BAM filtering
└── pNpS.py         # pN/pS analysis
resources/
└── genomad_db/     # geNomad database (not included)
```

## 3. Dependencies

This pipeline is built using **Snakemake**. All software dependencies are managed via **conda**.

* **Snakemake** (tested with v7.32.4)
* **geNomad** (tested with v1.8.1)
* **BWA** (tested with v0.7.17)
* **SAMtools** (tested with v1.19.2)
* **python** (tested with v3.11.9)
* **pandas** (tested with v2.2.0)
* **pysam** (tested with v0.22.1)
* **Biopython** (tested with v1.83)
* **LoFreq** (tested with v2.1.5)

## 4. System requirements

This pipeline is designed for high-performance computing (HPC) environments.

### Hardware requirements
* **Memory:** Minimum 32 GB RAM
* **CPU:** Scalable from 16 to 48+ cores

## 5. Installation & Setup

### 5.1 Clone the repository

```bash
git clone https://github.com/CSB5/SPMP_Phages.git
```

### 5.2 Install dependencies

```bash
conda create -n pnps -c conda-forge -c bioconda snakemake=7.32.4 genomad=1.8.1 bwa=0.7.17 samtools=1.19.2 pandas=2.2.0 pysam=0.22.1 biopython=1.83 lofreq=2.1.5
conda activate pnps
```

### 5.3 Download databases

* **[geNomad](https://portal.nersc.gov/genomad/quickstart.html#downloading-the-database):**
  The geNomad database path should be specified in `config.yaml`.

## 6. Usage

### 6.1 Inputs

The pipeline requires two input FASTA files. Specify their paths in `config.yaml`.

| File | Description |
| --- | --- |
| `v95.fna` | Input vOTUs |
| `decoy_mags.fna` | Decoy MAGs with prophages removed |

The pipeline also requires a set of paired-end Illumina sequencing read files as input. Modify `sample_id`, `R1_path` and `R2_path` in `samples.tsv` accordingly.

The following settings are defined in `config.yaml`:
* Reference FASTA (`reference_fasta`)
* Decoy FASTA (`decoy_fasta`)
* Output filename prefix (`output_prefix`)

### 6.2 Parameters

The following parameters are defined in `config.yaml`:
* Minimum % identity of read alignment (`min_id`, default: 95)
* Minimum % coverage of read (`min_cov`, default: 80)
* Whether to check if reads are properly paired (`check_proper_pair`, default: 1)
* Minimum % coverage breadth for vOTU detection in a sample (`covthresh`, default: 70)

### 6.3 Running

Navigate to the project subdirectory to run the pipeline:
```bash
cd pipelines/pNpS
snakemake -s pNpS.smk --cores 48
```

## 7. Outputs

All results are saved in the `results/` directory. The primary outputs are:

* **`<output_prefix>_vcf.tsv`:** Augmented variant calling TSV file, with columns:

  | sample_id | CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | AF | gene | codon_starting_index | snp_position_in_codon | ref_codon | mut_codon | mutation_type |
  | --------- | ----- | --- | -- | --- | --- | ---- | ------ | ---- | -- | ---- | -------------------- | --------------------- | --------- | --------- | ------------- |
  | sample ID | reference (vOTU) ID | SNP coordinate | `.` | reference allele | alternate allele | Phred quality score | `PASS` | VCF's `INFO` fields | allele frequency | gene ID | codon starting index (coordinate - 1) | SNP coordinate in codon | reference codon | mutant codon | `Intergenic`/`Syn`/`Non-syn`/`Stop` |

* **`<output_prefix>_gene_pNpS.tsv`:** Augmented gene annotation TSV file for vOTUs detected in each sample, with columns:
  
  | sample_id | CHROM | ... | possible_S | possible_NS | num_SNP_pos | S | NS | S_pseudocount | NS_pseudocount | pNpS |
  | --------- | ----- | --- | ---------- | ----------- | ----------- | - | -- | ------------- | -------------- | ---- |
  | sample ID | reference (vOTU) ID | geNomad's [`genes.tsv`](https://portal.nersc.gov/genomad/quickstart.html#understanding-the-outputs) (20 columns) | number of synonymous single-nucleotide substitutions | number of nonsynonymous single-nucleotide substitutions | number of SNP positions in gene+sample | number of synonymous SNPs | number of nonsynonymous SNPs | number of synonymous SNPs (+1 pseudocount) | number of nonsynonymous SNPs (+1 pseudocount) | pN/pS |

## 8. Citation

If you use this pipeline in your research, please cite:

* Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711
* **geNomad**: Camargo, A. P. et al. Identification of mobile genetic elements with geNomad. *Nature Biotechnology* **42**, 1303–1312 (2024).
* **BWA**: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics* **25**, 1754–1760 (2009).
* **SAMtools**: Danecek, P. et al. Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
* **LoFreq**: Wilm, A. et al. LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. *Nucleic Acids Research* **40**, 11189–11201 (2012).
