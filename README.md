# TWAS_ASSOC
>Implement the eQTL analysis (first step of TWAS) in Python.

![LOGO](/Image/TWASFUSION.png)
## Introduction
  `TWAS_ASSOC` is developed to employ FUSION, an advanced suite of tools for performing large-scale Transcriptome-Wide and Regulome-Wide Association Studies (TWAS and RWAS). 
 
  This implementation leverages the innovative work by [Gusev et al. (2016)](https://www.nature.com/articles/ng.3506) that integrates gene expression measurements with summary association statistics from large-scale GWAS. 
  
  The objective is to integrate genotype, gene expression and phenotype to gain insights into the genetic basis of complex traits

![Project Scheme](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fng.3506/MediaObjects/41588_2016_Article_BFng3506_Fig1_HTML.jpg?as=webp)

## Steps
- Preprocess: Prioritizing GTEx genes with h2 p-value < 0.05, particularly in Whole Blood relevance within the 1Mb gene region, following regression analysis for covariate removal.
  ```

  ```
- Step 1: Predicting gene expression utilizing cis-SNPs within a 1Mb range as predictors, employing single-gene regression for effect size estimation and employing FUSION methods.
  ```

  ```
- Step 2: Use estimated SNP weights to construct predicted expression and test association with GWAS traits.
  ```

  ```
  

## Installation
`TWAS_ASSOC` is available on , you can install the latest version via ``:
```commandline

```
