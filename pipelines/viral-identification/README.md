# Viral Identification Pipeline

This is a bioinformatics pipeline for identifying non-redundant viral populations from metagenomic assemblies leveraging **VirSorter2**, **CheckV**, and **geNomad**. Predicted viral sequences are filtered, aggregated across samples, and dereplicated at 95% ANI over 85% aligned fraction to generate a set of viral operational taxonomic units (vOTUs) for downstream analyses.

---

## 1. Overview

This pipeline performs the following steps:

1. **Viral prediction (VirSorter2)**
2. **Quality assessment and host trimming (CheckV):** Further host trimming (following [this SOP](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=3)).
3. **Filtering (VirSorter2+CheckV):** Viral sequences above a minimum length (5kbp) and viral ≥ host genes are retained.
4. **Viral prediction (geNomad)**
5. **Filtering (geNomad):** Viral sequences above a minimum length (5kbp) and below an FDR threshold (1%) are retained.
6. **Aggregation:** Viral sequences are combined across all samples.
7. **Dereplication:** Viral sequences are dereplicated at 95% ANI over 85% aligned fraction to generate a set of non-redundant vOTUs.

---

## 2. Directory structure

```text
viral-identification/
├── README.md
├── vir-id.smk
├── config/
│   └── config.yaml
├── data/
│   └── <sample.fna>  # metagenomic assemblies (not included)
├── resources/
│   └── genomad_db/   # geNomad database (not included)
└── results/          # outputs (not included)
```

---

## 3. Dependencies

This pipeline is built using **Snakemake**. All software dependencies are managed via **conda**.

* **Snakemake** (tested with v5.26.0)
* **VirSorter2** (tested with v2.2.4)
* **CheckV** (tested with v1.0.1)
* **geNomad** (tested with v1.8.1)
* **CD-HIT** (tested with v4.8.1)
* **python** (tested with v3.10.0)
* **pandas** (tested with v2.2.0)
* **Biopython** (tested with v1.83)

---

## 4. System requirements

This pipeline is designed for high-performance computing (HPC) environments.

### Hardware requirements
* **Memory:** Minimum 64 GB RAM
* **CPU:** Scalable from 16 to 48+ cores

---

## 5. Installation & Setup

### 5.1 Clone the repository

```bash
git clone https://github.com/CSB5/SPMP_Phages.git
cd viral-identification
```

### 5.2 Install dependencies

Create and activate a new conda environment:

```bash
conda create -n vir-id -c conda-forge -c bioconda snakemake=5.26.0 virsorter=2.2.4 checkv=1.0.1 genomad=1.8.1 cd-hit=4.8.1 pandas=2.2.0 biopython=1.83
conda activate vir-id
```

### 5.3 Download databases

* **VirSorter2 database:**
  See [https://github.com/jiarong/VirSorter2#download-database-and-dependencies](https://github.com/jiarong/VirSorter2#download-database-and-dependencies). Make sure to initialize the database using `virsorter setup` or `virsorter config`.

* **CheckV database:**
  See [https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database](https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database).

* **geNomad database:**
  See [https://portal.nersc.gov/genomad/quickstart.html#downloading-the-database](https://portal.nersc.gov/genomad/quickstart.html#downloading-the-database). The geNomad database path should be specified in `config/config.yaml`.

---

## 6. Usage

### 6.1 Inputs

The pipeline requires a set of metagenomic assembly files as input. The following settings are defined in `config/config.yaml`:

* Number of samples (`num_samples`)
* Sample name prefix (`sample_prefix`)
* Sample number padding (`padding`)
* Sample directory (`contigs_dir`)
* Sample filename suffix (`contigs_filename`)
* Contig header prefix (`contig_name_prefix`)

Currently, sample names are generated as:

```
<sample_prefix>1, <sample_prefix>2, ..., <sample_prefix><num_samples>
```

Modify parameters in `config.yaml` and/or `SAMPLES` in `vir-id.smk` for other sample naming conventions.

### 6.2 Parameters

The following parameters are defined in `config/config.yaml`:

* Minimum viral sequence length (`min_len`, default: 5000bp)
* VirSorter2 minimum score (`min_vs2_score`, default: 0.5)
* geNomad maximum FDR (`max_gmd_fdr`, default: 0.01)
* geNomad database path (`genomad_db`)

### 6.3 Running

To execute the pipeline, run the following command from the project root:

```bash
snakemake -s vir-id.smk --cores 48
```

---

## 7. Outputs

All results are saved in the `results/` directory. The primary outputs are:

* **`viral.fna`**
  Viral sequences identified by VirSorter2+CheckV and geNomad, aggregated across all samples, and dereplicated at 100% ANI.

* **`v95.fna`**
  Non-redundant vOTUs at 95% ANI over 85% aligned fraction.

These outputs are used as inputs for downstream viral clustering, host association and other analyses.

---

## 8. Authors & Citation

Author: [Hanrong Chen](mailto:chenhr@a-star.edu.sg), Genome Institute of Singapore, A\*STAR

If you use this pipeline in your research, please cite:

> Chen, H. et al. [...] (2026).
See Supplementary Fig. 1 for a schematic overview of this pipeline.

This pipeline utilizes several third-party tools. Please also cite the following primary publications:

* **VirSorter2**: Guo, J. et al. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome* **9**, 37 (2021).
* **CheckV**: Nayfach, S. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology* **39**, 578–585 (2021).
* **geNomad**: Camargo, A. P. et al. Identification of mobile genetic elements with geNomad. *Nature Biotechnology* **42**, 1303–1312 (2024).
* **CD-HIT**: Fu, L. et al. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* **28**, 3150–3152 (2012).
