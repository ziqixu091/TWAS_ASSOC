# TWAS_ASSOC

## Introduction
  `TWAS_ASSOC` is developed to employ FUSION TWAS, an advanced suite of tools for performing large-scale Transcriptome-Wide Association Studies created by [Gusev et al. Nature Genetics 2016](https://www.nature.com/articles/ng.3506). In this project, we aim to build the gene model and conduct eQTL (expression quantitative trait loci) mapping, the first step of TWAS, using publicly available genotype and LCL gene expression dataset [Geuvadis](https://www.internationalgenome.org/data-portal/data-collection/geuvadis) that contains 300 samples of the European population and 89 samples of the African population. This step aims to link the gene expressions with the genetic markers to identify significant gene-SNP pairs that could explain the variation of expression levels. Because of the long time span to get access to publicly available GWAS (Genome-wide Association Studies) genotype and phenotype data, we will only focus on the first half of the TWAS model. We replicate `lasso`, `Ridge regression`, `Elastic net`, `marginal z-scores`, and `Bayesian LMM` models for conducting the association tests and provide user-friendly software implemented in Python. 


## Installation
To set up the appropriate environment:
```

```
To install the software, simply do the following:
```
git clone https://github.com/ziqixu091/TWAS_ASSOC.git
```
If you are using the BSLMM model to conduct association tests, download and install the [GEMMA](https://xiangzhou.github.io/software/) software, and add it to path. 
We implemented the other models using Python packages. Please refer to Analysis and Output for more details.

## Analysis and Outputs
- Data: The Geuvadis data is not provided in this repository; however, it could be easily downloaded on this [page](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1?query=GEUVADIS). The preprocessing of data involves using [PLINK](https://www.cog-genomics.org/plink/) tools where the detailed process is described in the `preprocess.txt` file. 

- Input: For different models, there are different input formats
  * lasso:
  ```

  ```
  * Ridge regression:
  ```

  ```
  * Elastic Net:
  ```

  ```
  * Top 1 - marginal z-scores:
  ```

  ```
  * Bayesian LMM:
  ```

  ```
  
- Output: The software will output both the effect size of the linear regression and also the r^2 statistics.
  
- Analysis: We focus on 10,000 disease-important GWAS genes and for each gene, we look at the SNPs within the 1Mb cis-window. We also conduct cis-heritability estimation using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) for each gene for all the European samples. We will use the metrics of adjusted accuracy (r^2/cis-h2) to compare the performance between models. 


