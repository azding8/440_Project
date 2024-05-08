## Overview
This repo contains the necessary scripts and instructions to reproduce the figures in the 20.440 final written report.

Endometriosis is a common and morbid disease that affects approximately 10-15% of reproductive-age women worldwide. 
The disease is characterized by deposits of endometrium-like tissue outside of the endometrium lining which grows and sheds throughout the menstrual cycle.
Despite its prevalence, diagnosis remains challenging and can take up to 10 years since symptoms are often similar to other pathologies such as cancer. 
The current standard of cancer diagnosis is highly accurate where well-defined biomarkers from biopsies are used to inform specific treatment regimens for unique cancer subtypes (triple negative, HER2 positive, etc.). 
The current standard of endometriosis diagnosis is limited to laparoscopy since there are no well-defined biomarkers for endometriosis subtypes.
Our study aims to identify candidate diagnostic biomarkers of the most common form of endometriosis (endometrioma) to allow the incorporation of biopsies into the standard endometriosis diagnostic workflow.
In aim 1 of our project, we leverage publicly available single-cell RNA sequencing (scRNA-seq) data from Fonseca et al. [1] to characterize transcriptomic differences between endometrioma and healthy ovarian epithelial cells.

## Data
In this study, we analyze the transcriptome of ovarian epithelial cells using a single-cell transcriptomic atlas of human endometriosis generated by Fonseca et al. 2023 [1]. 
This dataset profiled the transcriptomics of >370,000 individual cells from endometriomas (n=8), endometriosis (n=28), unaffected ovary (n=4), and endometriosis-free peritoneum (n=4) from a total of 21 different patients to generate a cellular atlas of endometrial-type epithelial cells, stromal cells, and additional cell populations. 
We extend the analysis performed in Fonseca et al. [1] by delving deeper into the subset of epithelial cells in endometrioma and unaffected ovaries. 

Processed Seurat objects are available at the following link:
https://cedars.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02

Primary sequencing data are available at the following link:
[GSE213216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216)

## Folder structure
`Cluster_scripts` contains shell scripts used to download fastqs from Sequence Read Archive (SRA) and process them into spliced/unspliced count matrices.

`Figures` contains PNGs of the figures in the written report.

`R_scripts` contains the R script used to perform analysis on ovarian epithelial cells.

`CSVs` contains a CSV of the DEGs between endometrioma and unaffected ovary epithelial cells.

## Installation
How do I run your code?
1. Download the required data files and unzip them
4. Modify the paths in `Ovary_Endometrioma_Analysis.R` to match your preferred local file structure
5. Run the code (make sure to install libraries which you are missing)

See comments in Rscript for required libraries and versions.

## References

1. Fonseca MAS, Haro M, Wright KN, Lin X et al. Single-cell transcriptomic analysis of endometriosis. Nat Genet 2023 Feb;55(2):255-267. PMID: 36624343

