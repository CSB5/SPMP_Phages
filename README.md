# SPMP Phages

Scripts used to generate the results and figures for the manuscript titled: "GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome".

**Project status (29 Jan 2026):** This repository is being actively curated for manuscript submission.
* **Available:** Core computational pipelines used to generate key results of the paper (viral OTUs, viral family-level clusters, prevalence and abundance estimation).
* **In progress:** Additional analysis pipelines and figure-generation notebooks, to be released shortly.

## Available pipelines

* [Identifying viral OTUs from metagenomic assemblies](./pipelines/viral-identification/README.md)
* [Constructing viral family-level clusters](./pipelines/vfc-construction/README.md)
* [Estimating prevalence and abundance in metagenomic samples](./pipelines/coverage-estimation/README.md)

## Raw data availability

* All Supplementary Data files, identified viral sequences, and representative vOTU sequences are available via Zenodo at [DOI:10.5281/zenodo.18253940](https://zenodo.org/records/18253940).
* Hi-C (n=84) and VLP (n=64) metagenomic sequencing reads are available from the European Nucleotide Archive (ENA) under project accession [PRJEB106095](https://www.ebi.ac.uk/ena/browser/view/PRJEB106095). Illumina (n=109), Oxford Nanopore (n=109), and Hi-C (n=24) metagenomic sequencing reads are available under project accession [PRJEB49168](https://www.ebi.ac.uk/ena/browser/view/PRJEB49168).

## Contact

For questions regarding the code and data in this repository, please contact:
[Hanrong Chen](mailto:chenhr@a-star.edu.sg), [Niranjan Nagarajan](mailto:niranjan@nus.edu.sg)

## Citation

Chen, H. et al. GuFi phages represent the most prevalent viral family-level clusters in the human gut microbiome. bioRxiv. https://doi.org/10.64898/2026.01.26.701711
