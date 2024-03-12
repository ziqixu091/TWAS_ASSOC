import numpy as np
from helper_fns import *
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

from multiprocessing import Pool

def weights_marginal(genos, pheno, beta=False):
    """
    Calculates the marginal effect weights based on genotype (genos) and phenotype (pheno) data.
    If beta is True, adjusts the calculation by dividing by (n - 1),
    otherwise divides by the square root of (n - 1) for scaling.

    Parameters:
    - genos: 2D numpy array of genotype data (individuals x genetic variants)
    - pheno: 1D numpy array of phenotype data
    - beta: Boolean, if True adjusts the calculation differently

    Returns:
    - eff_wgt: 1D numpy array of effect weights
    """
    n = genos.shape[0]
    if beta:
        eff_wgt = np.sum(genos * pheno[:, np.newaxis], axis=0) / (n - 1)
    else:
        eff_wgt = np.sum(genos * pheno[:, np.newaxis], axis=0) / np.sqrt(n - 1)
    return eff_wgt


def Marginal(X, Y, beta=False, test_size=0.2, random_state=42):
    X_std = StandardScaler().fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X_std, Y, test_size=test_size, random_state=random_state)
    eff_wgt = weights_marginal(X_train, Y_train, beta)
    eff_wgt_squared = eff_wgt ** 2
    max_idx = np.argmax(eff_wgt_squared)
    # Keep only the one with the highest squared effect weight and set the rest to 0
    eff_wgt_filtered = np.zeros_like(eff_wgt)
    eff_wgt_filtered[max_idx] = eff_wgt[max_idx]
    return {
        "eff_wgt": eff_wgt,
        "eff_wgt_filtered": eff_wgt_filtered,
        # "X_train": X_train,
        # "X_test": X_test,
        # "Y_train": Y_train,
        # "Y_train_pred": np.dot(X_train, eff_wgt),
        # "Y_test": Y_test,
        # "Y_test_pred": np.dot(X_test, eff_wgt),
        "r2_train": r2_score(Y_train, np.dot(X_train, eff_wgt_filtered)),
        "r2_test": r2_score(Y_test, np.dot(X_test, eff_wgt_filtered))
    }


def process_gene(gene_id, protein_genes, ancsetry, y_full_df, test_size, bfile, random_state):
    try:
        processed_geno, X, Y = process_one_gene(gene_id, protein_genes, ancsetry, y_full_df, bfile)
        snps = processed_geno["snp_info"]
        return gene_id, Marginal(X, Y, test_size=test_size, random_state=random_state), snps
    except ValueError:
        print("No snps for gene ", gene_id)
        return None

def Marginal_all_genes(protein_genes, y_full_df, ancsetry, bfile=None, test_size=0.2, random_state=42, save_snps=True):
    results = {}
    snps_list = []
    gene_ids = protein_genes["gene_id"]
    pool = Pool()
    results_list = pool.starmap(process_gene, [(gene_id, protein_genes, ancsetry, y_full_df, test_size, bfile, random_state) for gene_id in gene_ids])
    pool.close()
    pool.join()
    for result in results_list:
        if result is not None:
            gene_id, result_data, snps = result
            snps["effect_weights"] = result_data["eff_wgt"]
            snps_list.append(snps)
            if save_snps:
                os.makedirs(f"./project_data/results/effect/{ancsetry}", exist_ok=True)
                snps.to_csv(f"./project_data/results/effect/{ancsetry}/Marginal_{gene_id}_effect_weights.csv")
            results[gene_id] = result_data
    
    train_r2 = np.mean([results[gene_id]["r2_train"] for gene_id in results])
    test_r2 = np.mean([results[gene_id]["r2_test"] for gene_id in results])
    print(f"Average train R^2: {train_r2}, Average test R^2: {test_r2}")
    
    return results, snps_list


if __name__ == "__main__":
    EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes = load_data()

    if not os.path.exists("./project_data/results"):
        os.makedirs("./project_data/results")

    EUR_results = Marginal_all_genes(EUR_protein_genes, EUR_ge_regressed, "EUR")
    with open("./project_data/results/marginal_results_EUR.pkl", "wb") as f:
        pickle.dump(EUR_results, f)
    print("Results saved in ./project_data/results/marginal_results_EUR.pkl")

    YRI_results = Marginal_all_genes(YRI_protein_genes, YRI_ge_regressed, "YRI")
    with open("./project_data/results/marginal_results_YRI.pkl", "wb") as f:
        pickle.dump(YRI_results, f)
    print("Results saved in ./project_data/results/marginal_results_YRI.pkl")
