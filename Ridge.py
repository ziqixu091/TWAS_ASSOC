from helper_fns import *
from tqdm import tqdm
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

def Ridge(X, Y, alpha=0.1, test_size=0.2, random_state=42):
    X_std = StandardScaler().fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X_std, Y, test_size=test_size, random_state=random_state)
    clf = linear_model.Ridge(alpha=alpha)
    clf.fit(X_train, Y_train)
    return {
        "clf": clf,
        # "X_train": X_train,
        # "X_test": X_test,
        # "Y_train": Y_train,
        # "Y_train_pred": clf.predict(X_train),
        # "Y_test": Y_test,
        # "Y_test_pred": clf.predict(X_test),
        "r2_train": clf.score(X_train, Y_train),
        "r2_test": clf.score(X_test, Y_test)
    }

def Ridge_all_genes(protein_genes, y_full_df, ancsetry, bfile, alpha=0.1, test_size=0.2, random_state=42, save_snps=True):
    results = {}
    snps_list = []
    for gene_id in tqdm(protein_genes["gene_id"]):
        try:
            processed_geno, X, Y = process_one_gene(gene_id, protein_genes, ancsetry, y_full_df, bfile)
        except ValueError:
            print("No snps for gene ", gene_id)
            continue
        results[gene_id] = Ridge(X, Y, alpha=alpha, test_size=test_size, random_state=random_state)
        snps = processed_geno["snp_info"]
        snps["effect_weight"] = results[gene_id]["clf"].coef_
        snps_list.append(snps)
        if save_snps:
            os.makedirs(f"./project_data/results/effect/{ancsetry}", exist_ok=True)
            snps.to_csv(f"./project_data/results/effect/{ancsetry}/Ridge_{gene_id}_effect_weights.csv")
    
    train_r2 = np.mean([results[gene_id]["r2_train"] for gene_id in results])
    test_r2 = np.mean([results[gene_id]["r2_test"] for gene_id in results])
    print(f"Average train R^2: {train_r2}, Average test R^2: {test_r2}")
    
    return results, snps_list


if __name__ == "__main__":
    EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes = load_data()

    if not os.path.exists("./project_data/results"):
        os.makedirs("./project_data/results")

    EUR_results = Ridge_all_genes(EUR_protein_genes, EUR_ge_regressed, "EUR")
    with open("./project_data/results/ridge_results.pkl", "wb") as f:
        pickle.dump(EUR_results, f)
    print("Results saved in ./project_data/results/ridge_results.pkl")

    YRI_results = Ridge_all_genes(YRI_protein_genes, YRI_ge_regressed, "YRI")
    with open("./project_data/results/ridge_results.pkl", "wb") as f:
        pickle.dump(YRI_results, f)
    print("Results saved in ./project_data/results/ridge_results.pkl")