import subprocess
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import concurrent
from tqdm import tqdm
import os
from pyplink import PyPlink
from sklearn.metrics import r2_score

def weights_bslmm(input, bv_type, snp, out=None, gemma_path="gemma", sys_print=False):
    """
    Run BSLMM analysis using GEMMA and extract effect sizes for specified SNPs.

    Parameters:
    - input: Base name for input files.
    - bv_type: Specifies the type of BSLMM analysis.
    - snp: List or array of SNP identifiers for which weights are to be calculated.
    - out: Optional. Specifies the base name for output files. Defaults to None.
    - gemma_path: Path to the GEMMA executable. Defaults to 'gemma'.
    - sys_print: If True, prints the GEMMA command output.

    Returns:
    - A numpy array of effect weights for the input SNPs.
    """
    if out is None:
        out = f"{input}.BSLMM"

    # # Constructing the GEMMA command
    arg = f"{gemma_path} -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile {input} -bslmm {bv_type} -o {out}"

    # Execute the GEMMA command
    result = subprocess.run(arg, shell=True, capture_output=not sys_print)
    if not sys_print:
        print(result.stdout.decode())  # Optional: print GEMMA output for debugging.

    # Read the output parameter file
    try:
        eff = pd.read_table(f"./output/{out}.param.txt", header=0, sep='\t')
    except FileNotFoundError:
        raise FileNotFoundError("GEMMA output file not found. Check GEMMA execution and output path.")

    # Initialize effect weights with NaN for all SNPs
    eff_wgt = pd.Series(np.nan, index=snp)

    # Match SNPs and assign weights
    for i, snp_id in enumerate(snp):
        if snp_id in eff['rs'].values:
            row = eff.loc[eff['rs'] == snp_id].iloc[0]
            eff_wgt.at[snp_id] = row['alpha'] + row['beta'] * row['gamma']

    return eff_wgt


def train_pred_bslmm(gene_id, ancestry, bv_type=1, out=None, gemma_path="gemma", sys_print=False):
    os.chdir(f"/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene_id}")
    all_files = os.listdir()
    # Pick the one ends with _train.bed
    bfile_train = [file for file in all_files if file.endswith("_train.bed")][0].split(".")[0]
    bfile_test = [file for file in all_files if file.endswith("_test.bed")][0].split(".")[0]
    with PyPlink(bfile_train) as bed:
        # Get SNP names
        bim = bed.get_bim()
        snps = bim.index.tolist()
        fam = bed.get_fam()
        # Get effect weights
        weights = weights_bslmm(bfile_train, bv_type, snps, out, gemma_path, sys_print)
        weights_unfiltered = weights
        train_y = fam[["status"]].values.flatten()
        X = np.zeros((len(train_y), len(snps)))

        for snp_id, genotypes in bed.iter_geno_marker(snps):
            X[:, snps.index(snp_id)] = genotypes

        # Find index of nan in weights
        nan_idx = np.argwhere(np.isnan(weights)).flatten()
        # Remove nan from weights and X
        weights = np.delete(weights, nan_idx)
        X = np.delete(X, nan_idx, axis=1)
        train_pred_y = np.dot(X, weights)
        train_r2 = r2_score(train_y, train_pred_y)

    with PyPlink(bfile_test) as bed:
        # Get SNP names
        bim = bed.get_bim()
        snps = bim.index.tolist()
        fam = bed.get_fam()

        test_y = fam[["status"]].values.flatten()
        X = np.zeros((len(test_y), len(snps)))

        for snp_id, genotypes in bed.iter_geno_marker(snps):
            X[:, snps.index(snp_id)] = genotypes
        X = np.delete(X, nan_idx, axis=1)
        test_pred_y = np.dot(X, weights)
        test_r2 = r2_score(test_y, test_pred_y)

    return {
        "eff_weights": weights,
        "eff_weights_unfiltered": weights_unfiltered,
        "train_r2": train_r2,
        "test_r2": test_r2
    }


if __name__ == "__main__":
    # Use ProcessPoolExecutor to execute the tasks in parallel
    with ProcessPoolExecutor() as executor:
        # Using list comprehension to create and immediately start all tasks
        # tqdm is used to display progress
        all_genes_EUR = os.listdir("/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/EUR/")
        futures = [executor.submit(train_pred_bslmm, gene_id, "EUR") for gene_id in all_genes_EUR]
        # Ensure all tasks are complete before moving on
        results = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            results.append(future.result())
        with open("./project_data/results/BSLMM_results_EUR.pkl", "wb") as f:
            pickle.dump(results, f)
        
    with ProcessPoolExecutor() as executor:
        # Using list comprehension to create and immediately start all tasks
        # tqdm is used to display progress
        all_genes_YRI = os.listdir("/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/YRI/")
        futures = [executor.submit(train_pred_bslmm, gene_id, "YRI") for gene_id in all_genes_YRI]
        # Ensure all tasks are complete before moving on
        results = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            results.append(future.result())
        with open("./project_data/results/BSLMM_results_YRI.pkl", "wb") as f:
            pickle.dump(results, f)
