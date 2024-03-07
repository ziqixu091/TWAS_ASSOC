import sys
import os
import pickle
import json

from helper_fns import *

from ElasticNet import ElasticNet_all_genes
from LASSO import LASSO_all_genes
from Ridge import Ridge_all_genes
from Marginal import Marginal_all_genes

def save_results(results, ancsetry, model):
    if not os.path.exists("./project_data/results"):
        os.makedirs("./project_data/results")

    with open(f"./project_data/results/{model}_results_{ancsetry}.pkl", "wb") as f:
        pickle.dump(results, f)
    print(f"Results saved in ./project_data/results/{model}_results_{ancsetry}.pkl")


if __name__ == "__main__":
    with open("config.json", "r") as f:
        config = json.load(f)

    bfile_path = config["bfile_path"]
    gene_exp_regressed_path = config["gene_exp_regressed_path"]
    protein_gene_path = config["protein_gene_path"]
    model = config["model"]
    ancestry = config["ancestry"]
    num_genes = config["num_genes"]
    save_model = config["save_model"]
    # plink_path = config["plink_path"]
    # gemma_path = config["gemma_path"]

    ge_regressed_chr, protein_genes_chr = load_custom_data(protein_gene_path=protein_gene_path,
                                                            ge_regressed_path=gene_exp_regressed_path,
                                                            chr_num=1, num_genes=num_genes)
    
    if "ElasticNet" in model:
        print("\nRunning ElasticNet on {} genes and {} samples".format(protein_genes_chr.shape[0], ge_regressed_chr.shape[0]))
        results, snps_list = ElasticNet_all_genes(protein_genes_chr, ge_regressed_chr, ancestry, bfile=bfile_path, save_snps=save_model)
        if save_model:
            save_results(results, ancestry, "ElasticNet")
    if "LASSO" in model:
        print("\nRunning LASSO on {} genes and {} samples".format(protein_genes_chr.shape[0], ge_regressed_chr.shape[0]))
        results, snps_list = LASSO_all_genes(protein_genes_chr, ge_regressed_chr, ancestry, bfile=bfile_path, save_snps=save_model)
        if save_model:
            save_results(results, ancestry, "LASSO")
    if "Ridge" in model:
        print("\nRunning Ridge on {} genes and {} samples".format(protein_genes_chr.shape[0], ge_regressed_chr.shape[0]))
        results, snps_list = Ridge_all_genes(protein_genes_chr, ge_regressed_chr, ancestry, bfile=bfile_path, save_snps=save_model)
        if save_model:
            save_results(results, ancestry, "Ridge")
    if "Marginal" in model:
        print("\nRunning Marginal on {} genes and {} samples".format(protein_genes_chr.shape[0], ge_regressed_chr.shape[0]))
        results, snps_list = Marginal_all_genes(protein_genes_chr, ge_regressed_chr, ancestry, bfile=bfile_path, save_snps=save_model)
        if save_model:
            save_results(results, ancestry, "Marginal")
    
