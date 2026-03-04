# SPMP Phages

This repository contains core computational pipelines and analysis scripts used to generate the results and figures for the manuscript:
> Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711

## Core computational pipelines

* [Identifying viral OTUs from metagenomic assemblies](./pipelines/viral-identification/README.md)
* [Constructing viral family-level clusters](./pipelines/vfc-construction/README.md)
* [Estimating prevalence and abundance in metagenomic samples](./pipelines/coverage-estimation/README.md)
* [Predicting virus-host associations from Hi-C data](./pipelines/HiC-host-association/README.md)
* [Variant calling and pN/pS analysis](./pipelines/pNpS/README.md)

## Figure-generation scripts

* [notebooks](./figures/notebooks/): Python notebooks and R scripts to reproduce all figures in the manuscript (tested with python v3.10.11 and R v4.5.0). To reproduce the figures, clone this repository and execute the code from within this directory.

* [data](./figures/data/): Input data files. Before executing the code, also download the Supplementary Data files from Zenodo and save them in this directory.

## Raw data availability

* All Supplementary Data files, identified viral sequences, and representative vOTU sequences are available via Zenodo at [DOI:10.5281/zenodo.18780256](https://zenodo.org/records/18780256).
* The hybrid MAGs used for host association and read coverage analysis are available from the European Nucleotide Archive (ENA) under project accession [PRJEB49168](https://www.ebi.ac.uk/ena/browser/view/PRJEB49168?show=analyses).
* Hi-C (n=84) and VLP (n=64) metagenomic sequencing reads are available from the European Nucleotide Archive (ENA) under project accession [PRJEB106095](https://www.ebi.ac.uk/ena/browser/view/PRJEB106095). Illumina (n=109), Oxford Nanopore (n=109), and Hi-C (n=24) metagenomic sequencing reads are available under project accession [PRJEB49168](https://www.ebi.ac.uk/ena/browser/view/PRJEB49168?show=reads).

## Contact

For questions regarding the code and data in this repository, please contact [Hanrong Chen](mailto:chenhr@a-star.edu.sg) and [Niranjan Nagarajan](mailto:niranjan@nus.edu.sg).

## Citation

If you use this repository, please cite:
> Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711

Please also cite the lab's previous paper on the SPMP cohort:
> Gounot, J.-S. et al. Genome-centric analysis of short and long read metagenomes reveals uncharacterized microbiome diversity in Southeast Asians. *Nature Communications*, **13**, 6044 (2022).
