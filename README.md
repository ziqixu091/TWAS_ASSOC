# TWAS_ASSOC

## Introduction
  `TWAS_ASSOC` is developed to employ FUSION TWAS, an advanced suite of tools for performing large-scale Transcriptome-Wide Association Studies created by [Gusev et al. Nature Genetics 2016](https://www.nature.com/articles/ng.3506). In this project, we aim to build the gene model and conduct eQTL (expression quantitative trait loci) mapping, the first step of TWAS, using publicly available genotype and LCL gene expression dataset [Geuvadis](https://www.internationalgenome.org/data-portal/data-collection/geuvadis) that contains 300 samples of the European population and 89 samples of the African population. This step aims to link the gene expressions with the genetic markers to identify significant gene-SNP pairs that could explain the variation of expression levels. Because of the long time span to get access to publicly available GWAS (Genome-wide Association Studies) genotype and phenotype data, we will only focus on the first half of the TWAS model. We replicate `lasso`, `Ridge regression`, `Elastic net`, `marginal z-scores`, and `Bayesian LMM` models for conducting the association tests and provide user-friendly software implemented in Python. 


## Installation
To install the software, simply do the following:
```
git clone https://github.com/ziqixu091/TWAS_ASSOC.git
```
To set up the appropriate environment:
```
cd TWAS_ASSOC
```
```
conda env create -f env.yml
```
```
conda activate TWAS_ASSOC_284
```
If you are using the BSLMM model to conduct association tests, download and install the [GEMMA](https://xiangzhou.github.io/software/) software, and add it to path. 
We implemented the other models using Python packages. Please refer to Analysis and Output for more details.

## Demo Data
The demo data is hosted on Google Drive. To download it, please run the following command:
```
gdown --folder https://drive.google.com/drive/folders/1OS33asqrhRLkL3QM2qLHHzCN3G30cbaG?usp=drive_link
```

## Analysis and Outputs
- Data: The Geuvadis data is not provided in this repository; however, it could be easily downloaded on this [page](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1?query=GEUVADIS). The preprocessing of data involves using [PLINK](https://www.cog-genomics.org/plink/) tools where the detailed process is described in the `preprocess.txt` file. 

- Input: For different models and inputs, the `config.json` file needs to be modified accordingly. We've provided the `config.json` file, which runs all four models on the demo data. The `config.json` file contains the following parameters:

| Parameter                   | Default                                                     | Description                                           |
|-----------------------------|-------------------------------------------------------------|-------------------------------------------------------|
| **bfile_path**              | `./project_data/geno/YRI/GEUVADIS_YRI_chr1`                 | Path to the PLINK binary file                          |
| **protein_gene_path**       | `./project_data/GEUVADIS_YRI_protein_genes.tsv.gz`          | Path to the protein-coding gene information file      |
| **gene_exp_regressed_path** | `./project_data/GEUVADIS_YRI_ge_regressed.tsv.gz`           | Path to the gene expression file (covariates regressed) |
| **model**                   | `["ElasticNet", "LASSO", "Ridge", "Marginal"]`              | List of models to run                                 |
| **ancestry**                | `YRI`                                                       | Ancestry of the samples                               |
| **num_genes**               | `10`                                                        | Number of genes to run the model on                   |
| **save_model**              | `true`                                                      | Save the results of the model or not                  |

  
- Output: The software will output both the effect size of the linear regression and also the r^2 statistics in a pickle file. 
  
- Analysis: We focus on 10,000 disease-important GWAS genes and for each gene, we look at the SNPs within the 1Mb cis-window. We also conduct cis-heritability estimation using GCTA for each gene for all the European samples. We will use the metrics of adjusted accuracy (r^2/cis-h2) to compare the performance between models. 


