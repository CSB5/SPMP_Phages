# Predicting virus-host associations from Hi-C data

This is a bioinformatics pipeline for predicting virus-host associations from Hi-C sequencing data. Hi-C reads are mapped to an existing metagenomic assembly. The number of intra and inter-MAG linkages are compiled, and parameters for each MAG's contribution to Hi-C noise are inferred (see Methods of our [preprint](https://www.biorxiv.org/content/10.64898/2026.01.26.701711v1.full) for more details). P-values are computed for the observed number of virus-MAG linkages given the number of nonself-viral linkages and the inferred parameters for the MAG. An adjusted p-value < 0.05 was utilized for downstream analysis.

## 1. Overview

This pipeline performs the following steps:

1. **Read mapping:** Raw Hi-C reads are adapter-trimmed and quality-filtered using fastp, and mapped to an existing metagenomic assembly using BWA-MEM (options '-5SP').
2. **BAM filtering:** PCR duplicates are removed, and reads are filtered based on mapping identity (default: 95%) and coverage (default: 80%).
3. **Compilation of inter-MAG linkages:** The number of intra and inter-MAG Hi-C linkages are compiled.
4. **Solution of Hi-C noise model:** Parameters for each MAG's contribution to Hi-C noise are inferred.
5. **Compilation of virus-MAG linkages:** The number of virus-MAG Hi-C linkages are compiled. Care was taken to exclude viral self-linkages if the viral contig is binned in the MAG.
6. **P-value calculation:** P-values for the observed number of virus-MAG linkages are calculated given the total number of nonself-viral linkages and the inferred parameters for the MAG. P-values within each sample are adjusted for multiple testing using the Benjamini-Hochberg method.
7. **Aggregation of results:** Results are aggregated across samples.

## 2. Directory structure

```text
HiC-host-association/
├── README.md
├── hic.smk
├── config.yaml
├── samples.tsv                  # table of sample IDs and input path names
└── results/                     # outputs (not included)
data/
├── contigs/                     # metagenomic assemblies (not included)
├── HiC_reads/                   # directory containing Hi-C reads (not included)
├── mags/                        # directory containing TSV files with columns: contig_id, mag_id (not included)
└── viral/                       # directory containing TSV files with columns: vir_id, contig_id, start_bp, end_bp (not included)
scripts/
├── bamfilter.py                 # BAM filtering
├── bam_to_numlinkages.py        # compiles number of intra and inter-contig linkages
├── solve_hic_model.py           # solves Hi-C noise model for MAGs
└── bam_to_mapping_positions.py  # compiles linkage positions between contigs belonging to 2 specified lists
```

## 3. Dependencies

This pipeline is built using **Snakemake**. All software dependencies are managed via **conda**.

* **Snakemake** (tested with v7.32.4)
* **fastp** (tested with v0.23.4)
* **BWA** (tested with v0.7.17)
* **samblaster** (tested with v0.1.26)
* **SAMtools** (tested with v1.19.2)
* **python** (tested with v3.10.11)
* **pandas** (tested with v2.2.0)
* **numpy** (tested with v1.26.3)
* **pysam** (tested with v0.22.1)
* **scipy** (tested with v1.12.0)
* **statsmodels** (tested with v0.14.6)

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
conda create -n hic -c conda-forge -c bioconda snakemake=7.32.4 fastp=0.23.4 bwa=0.7.17 samblaster=0.1.26 samtools=1.19.2 python=3.10.11 pandas=2.2.0 numpy=1.26.3 pysam=0.22.1 scipy=1.12.0 statsmodels=0.14.6
conda activate hic
```

## 6. Usage

### 6.1 Inputs

The pipeline requires five input files for each sample. Specify their paths in `samples.tsv`.

| Column | Description |
| --- | --- |
| `sample_id` | Sample ID |
| `contigs_path` | Metagenomic assembly FASTA |
| `HiC_R1_path` | Hi-C R1 FASTQ |
| `HiC_R2_path` | Hi-C R2 FASTQ |
| `mag_contig_list` | TSV with columns: `contig_id`, `mag_id` |
| `viral_contig_list` | TSV with columns: `vir_id`, `contig_id`, `start_bp`, `end_bp` |

### 6.2 Parameters

The following parameters are defined in `config.yaml`:
* Minimum % identity of read alignment (`min_id`, default: 95)
* Minimum % coverage of read (`min_cov`, default: 80)
* Hi-C model hyperparameters (`hic_model`)

### 6.3 Running

Navigate to the project subdirectory to run the pipeline:
```bash
cd pipelines/HiC-host-association
snakemake -s hic.smk --cores 48
```

## 7. Outputs

All results are saved in the `results/` directory. The primary output is:

* **`virus-mag_HiC.tsv`:** TSV of virus-MAG pairs with at least one Hi-C linkage, with columns:

  | sample_id | vir_id1 | mag_id2 | num_linkages | Nself1 | M1 | N1 | Nself2 | M2 | max_othermag_m2 | N2 | e2 | p2 | pval | pval_corrected |
  |-----------|---------|---------|--------------|--------|----|----|--------|----|-----------------|----|----|----|------|----------------|
  | sample ID | viral sequence ID | MAG ID | no. of virus-MAG linkages | no. of viral self-linkages | no. of viral nonself-linkages | total no. of viral linkages | no. of MAG self-linkages | no. of MAG nonself-linkages | max no. of linkages between MAG and another MAG | total no. of MAG linkages | inferred `e` | inferred `p` | p-value | adjusted p-value |

## 8. Citation

If you use this pipeline in your research, please cite:

* Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711
* **fastp**:  Chen, S. fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta* **4**, e70078 (2025).
* **BWA**: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics* **25**, 1754–1760 (2009).
* **samblaster**: Faust, G. G. & Hall, I. M. SAMBLASTER: fast duplicate marking and structural variant read extraction. *Bioinformatics* **30**, 2503–2505 (2014).
* **SAMtools**: Danecek, P. et al. Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
