# SPMP Phages

This repository contains computational pipelines and analysis scripts used to generate the results and figures for the manuscript: "GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome".

## Available pipelines

Core computational pipelines used to generate the key results of the paper:

* [Identifying viral OTUs from metagenomic assemblies](./pipelines/viral-identification/README.md)
* [Constructing viral family-level clusters](./pipelines/vfc-construction/README.md)
* [Estimating prevalence and abundance in metagenomic samples](./pipelines/coverage-estimation/README.md)
* [Predicting virus-host associations from Hi-C data](./pipelines/HiC-host-association/README.md)
* [Variant calling and pN/pS analysis](./pipelines/pNpS/README.md)

## Figures

All code and input data to reproduce the figures in the manuscript.

* **`notebooks/`**
  Python notebooks and R scripts used to generate all figures.

  * Python notebooks were run with **Python 3.10.11**
  * R scripts were run with **R 4.5.0**
  
  To reproduce the figures, clone this repository and execute the code from within this directory.

* **`data/`**
  Input data files required by the scripts.
  In addition, download the Supplementary Data files from Zenodo ([https://zenodo.org/records/18253940](https://zenodo.org/records/18253940)) and save them in this directory before running the code.

## Raw data availability

* All Supplementary Data files, identified viral sequences, and representative vOTU sequences are available via Zenodo at [DOI:10.5281/zenodo.18253940](https://zenodo.org/records/18253940).
* The hybrid MAGs used for host association and read coverage analysis are available from the European Nucleotide Archive (ENA) under project accession [PRJEB49168](https://www.ebi.ac.uk/ena/browser/view/PRJEB49168?show=analyses).
* Hi-C (n=84) and VLP (n=64) metagenomic sequencing reads are available from the European Nucleotide Archive (ENA) under project accession [PRJEB106095](https://www.ebi.ac.uk/ena/browser/view/PRJEB106095). Illumina (n=109), Oxford Nanopore (n=109), and Hi-C (n=24) metagenomic sequencing reads are available under project accession [PRJEB49168](https://www.ebi.ac.uk/ena/browser/view/PRJEB49168?show=reads).

## Contact

For questions regarding the code and data in this repository, please contact:
[Hanrong Chen](mailto:chenhr@a-star.edu.sg), [Niranjan Nagarajan](mailto:niranjan@nus.edu.sg)

## Citation

Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711

Please also cite the lab's previous paper on the SPMP cohort:
Gounot, J.-S. et al. Genome-centric analysis of short and long read metagenomes reveals uncharacterized microbiome diversity in Southeast Asians. *Nature Communications*, **13**, 6044 (2022).
