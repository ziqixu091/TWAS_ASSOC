# TWAS_ASSOC

## Introduction
  `TWAS_ASSOC` is developed to employ FUSION TWAS, an advanced suite of tools for performing large-scale Transcriptome-Wide Association Studies created by [Gusev et al. (2016)](https://www.nature.com/articles/ng.3506). In this project, we aim to build the gene model and conduct eQTL (expression quantitative trait loci) mapping, the first step of TWAS, using publicly available genotype and LCL gene expression dataset [Geuvadis]() that contains 300 samples of the European population and 89 samples of the African population. This step aims to link the gene expressions with the genetic markers to identify significant gene-SNP pairs that could explain the variation of expression levels. Because of the long time span to get access to publicly available GWAS (Genome-wide Association Studies) genotype and phenotype data, we will only focus on the first half of the TWAS model. We replicate `lasso`, `Ridge regression`, `elastic net`, `marginal z-scores`, and `Bayesian LMM` models for conducting the association tests and provide user-friendly software implemented in Python. 


## Installation
`TWAS_ASSOC` is available on GitHub, you can install it by simply:
```
git clone https://github.com/ziqixu091/TWAS_ASSOC.git
```
If you are using the BSLMM model to conduct association tests, download and install the [GEMMA]() software, add to path. 

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
  
