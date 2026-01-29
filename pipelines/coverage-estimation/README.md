# Estimating prevalence and abundance in metagenomic samples

This is a bioinformatics pipeline for mean coverage and relative abundance estimation of vOTUs in a set of samples. It maps a sample's Illumina sequencing reads to input vOTUs and runs custom scripts for BAM filtering and clipped coverage computation. vOTUs passing a coverage breadth threshold (default: 70%) are considered present in a sample, and their relative abundances (based on mean coverage) are computed.

An identical pipeline was run on non-viral contigs, used for determining thresholds for viral replication in viral-enriched samples.

A similar pipeline was constructed to estimate MAG mean coverages, used for vOTU versus host MAG coverage ratios. Here,
* Only MAGs assembled from the same sample were used as reference.
* A coverage breadth threshold was not set.

## 1. Overview

This pipeline performs the following steps:

1. **Read mapping:** Paired-end Illumina reads are mapped to input vOTUs using BWA-MEM.
2. **BAM filtering:** Reads are filtered based on mapping identity, coverage, and proper pairing.
3. **Clipped coverage computation:** The mean coverage over each sequence is computed after clipping the top and bottom 10% of base coverage values.
4. **vOTU detection and abundance estimation:** vOTUs passing a coverage breadth threshold (default: 70%) are considered present, and their relative abundances are computed.
5. **Aggregation of results:** The mean coverage and relative abundance estimates are aggregated across samples.

## 2. Directory structure

```text
coverage-estimation/
├── README.md
├── cov-est.smk
├── config/
│   └── config.yaml
├── scripts/
│   ├── bamfilter.py             # BAM filtering
│   └── get_clipped_coverage.py  # clipped coverage computation
├── data/
│   ├── <v95.fna>                # input vOTUs (not included)
│   └── <reads/>                 # directory containing Illumina reads (not included)
└── results/                     # outputs (not included)
```

## 3. Dependencies

This pipeline is built using **Snakemake**. All software dependencies are managed via **conda**.

* **Snakemake** (tested with v7.32.4)
* **BWA** (tested with v0.7.17)
* **SAMtools** (tested with v1.19.2)
* **python** (tested with v3.10.13)
* **pandas** (tested with v2.2.0)
* **pysam** (tested with v0.22.1)

## 4. System requirements

This pipeline is designed for high-performance computing (HPC) environments.

### Hardware requirements

* **Memory:** Minimum 16 GB RAM
* **CPU:** Scalable from 8 to 24+ cores

## 5. Installation & Setup

### 5.1 Clone the repository

```bash
   git clone https://github.com/CSB5/SPMP_Phages.git
   cd pipelines/coverage-estimation
```

### 5.2 Install dependencies

```bash
conda create -n cov-est -c conda-forge -c bioconda snakemake=7.32.4 bwa=0.7.17 samtools=1.19.2 pandas=2.2.0 pysam=0.22.1
conda activate cov-est
```

## 6. Usage

### 6.1 Inputs

The pipeline requires a reference FASTA and paired-end Illumina sequencing read files as input.

The following settings are defined in `config/config.yaml`:

* Reference FASTA (`reference_fasta`)
* Prefix for id and len of references (`ref_prefix`)
* Final output filname prefix (`output_prefix`)
* Number of samples (`num_samples`)
* Sample name prefix (`sample_prefix`)
* Sample number padding (`padding`)
* Directory containing sequencing read files (`reads_dir`)
* R1 filename format (`reads_r1`)
* R2 filename format (`reads_r2`)

By default, sample names are generated as:
  ```
  <sample_prefix>1, <sample_prefix>2, ..., <sample_prefix><num_samples>
  ```

### 6.2 Parameters

The following parameters are defined in `config/config.yaml`:

* Minimum % identity of read alignment (`min_id`, default: 95)
* Minimum % coverage of read (`min_cov`, default: 80)
* Whether to check that reads are properly paired (`check_proper_pair`, default: 1)
* Minimum % coverage breadth for vOTU detection in a sample (`covthresh`, default: 70)

### 6.3 Running

To execute the pipeline, run the following command from the project root:

```bash
snakemake -s cov-est.smk --cores 24
```

## 7. Outputs

All results are saved in the `results/` directory. The primary output is:

* **`<output_prefix>_abundance.tsv`:**
  The coverage breadth, mean coverage, and relative abundance of detected vOTUs in all samples.

## 8. Citation

If you use this pipeline in your research, please cite:

* Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711
* **BWA**: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics* **25**, 1754–1760 (2009).
* **SAMtools**: Danecek, P. et al. Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
