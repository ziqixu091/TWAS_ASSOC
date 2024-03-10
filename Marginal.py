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
    eff_wgt_max = np.where(np.abs(eff_wgt) == np.max(np.abs(eff_wgt)), eff_wgt, 0)
    
    return {
        "eff_wgt": eff_wgt,
        "eff_wgt_max": eff_wgt_max,
        # "X_train": X_train,
        # "X_test": X_test,
        # "Y_train": Y_train,
        # "Y_train_pred": np.dot(X_train, eff_wgt),
        # "Y_test": Y_test,
        # "Y_test_pred": np.dot(X_test, eff_wgt),
        "r2_train": r2_score(Y_train, np.dot(X_train, eff_wgt_max)),
        "r2_test": r2_score(Y_test, np.dot(X_test, eff_wgt_max))
    }


def process_gene(gene_id, protein_genes, ancsetry, y_full_df, test_size, random_state):
    try:
        processed_geno, X, Y = process_one_gene(gene_id, protein_genes, ancsetry, y_full_df)
        return gene_id, Marginal(X, Y, test_size=test_size, random_state=random_state)
    except ValueError:
        print("No snps for gene ", gene_id)
        return None

def Marginal_all_genes(protein_genes, y_full_df, ancsetry, test_size=0.2, random_state=42):
    results = {}
    gene_ids = protein_genes["gene_id"]
    pool = Pool()
    results_list = pool.starmap(process_gene, [(gene_id, protein_genes, ancsetry, y_full_df, test_size, random_state) for gene_id in gene_ids])
    pool.close()
    pool.join()
    for result in results_list:
        if result is not None:
            gene_id, result_data = result
            results[gene_id] = result_data
    return results


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
