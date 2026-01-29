# Constructing viral family-level clusters

This is a bioinformatics pipeline for constructing viral clusters (VCs) via hierarchical clustering of a set of input vOTUs and reference vOTUs with known taxonomic assignments. It estimates viral genome completeness using CheckV, defines protein clusters (PCs), and performs UPGMA clustering on high-quality and complete vOTUs based on % shared distinct PCs.

Inspection of reference vOTU clustering at different distance thresholds led to viral family-level clusters (VFCs) being defined at 20% shared PCs. Medium-quality vOTUs are subsequently recruited into the VFCs if they share 30% PCs with members of that VFC, a threshold chosen to minimize ambiguous assignments to >1 VFC.

## 1. Overview

This pipeline performs the following steps:

1. **Completeness estimation:** Completeness estimation of input and reference vOTUs using CheckV.
2. **Gene identification:** Gene annotation of medium-quality and above (≥50% complete) vOTUs using geNomad.
3. **Protein clustering:** Construction of protein clusters (PCs) from protein sequences using MCL.
4. **Similarity calculation:** Computation of pairwise similarities (percentage of shared distinct PCs) between all pairs of medium-quality and above vOTUs.
5. **Hierarchical clustering:** UPGMA clustering of high-quality and above (≥90% complete) vOTUs only.
6. **Dendrogram cutting:** Cluster labels (VCs) are generated across a range of distance thresholds, which are manually inspected to define viral family-level clusters (VFCs).
7. **Recruitment:** Medium-quality (50-90% complete) vOTUs are recruited into VFCs if they share at least 30% PCs with members of that VFC.

## 2. Directory structure

```text
vfc-construction/
├── README.md
├── vfc.smk
├── scripts/
│   ├── get_pair_counts.py
│   └── perform_hierarchical_clustering.py
├── data/
│   ├── v95.fna         # input vOTU sequences (not included)
│   └── ref_fams95.fna  # reference vOTU sequences with known taxonomic assignments (not included)
├── resources/
│   └── genomad_db/     # geNomad database (not included)
└── results/            # outputs (not included)
```

## 3. Dependencies

This pipeline is built using **Snakemake**. All software dependencies are managed via **conda**.

* **Snakemake** (tested with v5.26.0)
* **CheckV** (tested with v1.0.1)
* **geNomad** (tested with v1.8.1)
* **seqtk** (tested with v1.4)
* **DIAMOND** (tested with v2.1.8)
* **MCL** (tested with v22.282)
* **python** (tested with v3.10.14)
* **pandas** (tested with v2.2.0)
* **numpy** (tested with v1.26.3)
* **scipy** (tested with v1.12.0)

## 4. System requirements

This pipeline is designed for high-performance computing (HPC) environments.

### Hardware requirements
* **Memory:** Minimum 32 GB RAM
* **CPU:** Scalable from 8 to 48+ cores

### Resource utilization
* **Dataset:** >30,000 input and >19,000 reference vOTUs
* **Hardware:** 48 CPUs on a high-compute node
* **Max memory (Peak RSS):** ~22 GB
* **Wall-clock time:** ~3 hours

## 5. Installation & Setup

### 5.1 Clone the repository

```bash
   git clone https://github.com/CSB5/SPMP_Phages.git
   cd pipelines/vfc-construction
```

### 5.2 Install dependencies

```bash
conda create -n vfc -c conda-forge -c bioconda snakemake=5.26.0 checkv=1.0.1 genomad=1.8.1 seqtk=1.4 diamond=2.1.8 mcl=22.282 pandas=2.2.0 numpy=1.26.3 scipy=1.12.0
conda activate vfc
```

### 5.3 Download databases

* **[CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database)**

* **[geNomad](https://portal.nersc.gov/genomad/quickstart.html#downloading-the-database):**
  Save `genomad_db/`, or a symlink to it, in the `resources/` directory.

## 6. Usage

The pipeline requires two primary input FASTA files. Save them in the `data/` directory.

| File | Description |
| --- | --- |
| `v95.fna` | Input vOTUs |
| `ref_fams95.fna` | Reference vOTUs with known taxonomic assignments |

To execute the pipeline, run the following command from the project root:

```bash
snakemake -s vfc.smk --cores 48
```

## 7. Outputs

All results are saved in the `results/` directory. The primary output is:

* **`vr95_c50-upgma.tsv`:**
  A tab-separated table containing viral cluster (VC) labels for all ≥90% complete and some 50-90% complete input and reference vOTUs, across a range of distance thresholds. Columns like `VC60` represent VCs formed at a 60% distance threshold (40% shared distinct PCs). These hierarchical clusters enable comparison with reference clade assignments, facilitating analysis of viral groups at the genus, family, and order levels.

## 8. Citation

If you use this pipeline in your research, please cite:

* Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711
* **CheckV**: Nayfach, S. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology* **39**, 578–585 (2021).
* **geNomad**: Camargo, A. P. et al. Identification of mobile genetic elements with geNomad. *Nature Biotechnology* **42**, 1303–1312 (2024).
* **DIAMOND**: Buchfink, B. et al. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods* **18**, 366–368 (2021).
* **MCL**: Van Dongen, S. Graph Clustering Via a Discrete Uncoupling Process. *SIAM J. Matrix Anal. Appl.* **30**, 121–141 (2008).
