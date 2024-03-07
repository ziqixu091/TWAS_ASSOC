import pandas as pd
import numpy as np

import os
import sys
from typing import List, Tuple, Dict, Any, Optional
import pickle

try:
    os.chdir('/Users/xbh0403/Desktop/TWAS_ASSOC')
except:
    os.chdir('/new-stg/home/banghua/TWAS_ASSOC')

from pyplink import PyPlink

def load_data():
    EUR_ge_regressed = pd.read_csv("./project_data/GEUVADIS_EUR_ge_regressed.tsv.gz", sep="\t", index_col=0, compression="gzip")
    YRI_ge_regressed = pd.read_csv("./project_data/GEUVADIS_YRI_ge_regressed.tsv.gz", sep="\t", index_col=0, compression="gzip")

    EUR_protein_genes = pd.read_csv("./project_data/GEUVADIS_EUR_protein_genes.tsv.gz", sep="\t", index_col=0, compression="gzip")
    EUR_protein_genes["chr"] = EUR_protein_genes.index
    EUR_protein_genes.reset_index(drop=True, inplace=True)
    YRI_protein_genes = pd.read_csv("./project_data/GEUVADIS_YRI_protein_genes.tsv.gz", sep="\t", index_col=0, compression="gzip")
    YRI_protein_genes["chr"] = YRI_protein_genes.index
    YRI_protein_genes.reset_index(drop=True, inplace=True)

    print("Shapes of the dataframes:", EUR_ge_regressed.shape, YRI_ge_regressed.shape, EUR_protein_genes.shape, YRI_protein_genes.shape)

    return EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes

def find_snps_in_gene(chr_num: int, start: int, end: int, ancsetry: str) -> pd.DataFrame:
    """
    Find SNPs within a specified genomic region.

    Parameters:
    chr_num (int): Chromosome number.
    start (int): Start position of the genomic region.
    end (int): End position of the genomic region.
    bim (pd.DataFrame): DataFrame containing SNP information.

    Returns:
    pd.DataFrame: DataFrame containing SNPs within the specified genomic region.
    
    """
    with PyPlink("./project_data/geno/"+ancsetry+"/GEUVADIS_"+ancsetry+"_chr"+str(chr_num)) as bed:
        # Getting the BIM and FAM
        bim = bed.get_bim()
        bim["snp"] = bim.index
        bim.reset_index(drop=True, inplace=True)
        bim = bim[["snp", "chrom", "pos", "cm", "a1", "a2"]]
    return bim[(bim["pos"] >= start) & (bim["pos"] <= end)]

def process_geno(ancsetry: str, chr_num: int, snps: List[str], start: int, gene_id: str, save_result: bool = False) -> np.ndarray:
    """
    Process the genotype data.

    Parameters:
    ancsetry (str): Ancestry of the individual.
    chr_num (int): Chromosome number.
    individual (str | None): Individual ID.
    snps (List[str]): List of SNPs.

    Returns:
    pd.DataFrame: Processed genotype data.
    """
    if len(snps) == 0:
        raise ValueError("No SNPs provided.")
    
    with PyPlink("./project_data/geno/"+ancsetry+"/GEUVADIS_"+ancsetry+"_chr"+str(chr_num)) as bed:
        # Getting the BIM and FAM
        bim = bed.get_bim()
        bim["snp"] = bim.index
        bim.reset_index(drop=True, inplace=True)
        bim = bim[["snp", "chrom", "pos", "cm", "a1", "a2"]]
        snp_info = bim[bim["snp"].isin(snps)]

        fam = bed.get_fam()
        iids = fam["iid"].tolist()

        results = {
            "snp_info": snp_info,
            "iids": iids,
        }

        keep_bool = [1] * len(iids)

        for snp_id, genotypes in bed.iter_geno_marker(snps):
            results[snp_id] = genotypes
            genotypes_kept = [0 if x == -1 else 1 for x in genotypes]
            keep_bool = [x*y for x, y in zip(keep_bool, genotypes_kept)]

        results["keep_bool"] = keep_bool

        for snp_id, genotypes in bed.iter_geno_marker(snps):
            results[snp_id+"_filtered"] = [x for x, y in zip(genotypes, keep_bool) if y == 1]
        
        results["iids_filtered"] = [x for x, y in zip(iids, keep_bool) if y == 1]

        X = np.zeros((len(results["iids_filtered"]), snp_info.shape[0]))

        for i in range(snp_info.shape[0]):
            X[:, i] = results[snps[i]+"_filtered"]

        if save_result:
            with open("./project_data/processed_Xy/"+gene_id+"_results.pkl", "wb") as f:
                pickle.dump(results, f)
            # save x as npy
            np.save("./project_data/processed_Xy/X/"+gene_id+"_X.npy", X)

        return results, X
    
def get_y(gene_id: str, iids: List[str], y_full_df: pd.DataFrame, save_y: bool = False) -> np.ndarray:
    """
    Get the Y values for a specific gene and ancestry.

    Parameters:
    gene_id (str): Gene ID.
    iids (List[str]): List of individual IDs.
    y_full_df (pd.DataFrame): DataFrame containing the Y values.

    Returns:
    np.ndarray: Array of Y values for the specified gene and ancestry.
    """
    all_columns = y_full_df.columns
    id_gene = all_columns.get_loc(gene_id)
    all_individuals = y_full_df.index
    iids_keep = [True if x in iids else False for x in all_individuals]
    y = y_full_df.iloc[iids_keep, id_gene].values
    if save_y:
        np.save("./project_data/processed_Xy/y/"+gene_id+"_y.npy", y)
    return y

def process_one_gene(gene_id: str, protein_genes: pd.DataFrame, ancsetry: str, y_full_df: pd.DataFrame) -> np.ndarray:
    gene = protein_genes[protein_genes["gene_id"] == gene_id]
    assert gene.shape[0] == 1
    chr_num = gene["chr"].values[0]
    start = gene["start"].values[0]
    end = gene["end"].values[0]
    gene_name = gene["name"].values[0]
    # print("Gene id: ", gene_id, " Gene name: ", gene_name, " Chr: ", chr_num, " Start: ", start, " End: ", end)
    
    snps = find_snps_in_gene(chr_num, start, end, ancsetry)
    snps_name = snps["snp"].tolist()
    # print("Number of SNPs: ", len(snps_name))
    processed_geno, X = process_geno(ancsetry, chr_num, snps_name, start, gene_id, save_result=False)
    y = get_y(gene_id, processed_geno["iids_filtered"], y_full_df, save_y=False)
    return processed_geno, X, y
