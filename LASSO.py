#%%
from helper_fns import *
from tqdm import tqdm
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

def LASSO(X, Y, alpha=0.1, test_size=0.2, random_state=42):
    X_std = StandardScaler().fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X_std, Y, test_size=test_size, random_state=random_state)
    clf = linear_model.Lasso(alpha=alpha)
    clf.fit(X_train, Y_train)
    return {
        "clf": clf,
        "X_train": X_train,
        "X_test": X_test,
        "Y_train": Y_train,
        "Y_train_pred": clf.predict(X_train),
        "Y_test": Y_test,
        "Y_test_pred": clf.predict(X_test),
        "r2_train": clf.score(X_train, Y_train),
        "r2_test": clf.score(X_test, Y_test)
    }

def LASSO_all_genes(protein_genes, y_full_df, ancsetry, alpha=0.1, test_size=0.2, random_state=42):
    results = {}
    for gene_id in tqdm(protein_genes["gene_id"]):
        try:
            processed_geno, X, Y = process_one_gene(gene_id, protein_genes, ancsetry, y_full_df)
        except ValueError:
            print("No snps for gene ", gene_id)
            continue
        results[gene_id] = LASSO(X, Y, alpha=alpha, test_size=test_size, random_state=random_state)
    return results

#%%
if __name__ == "__main__":
    EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes = load_data()

    if not os.path.exists("./project_data/results"):
        os.makedirs("./project_data/results")

    EUR_results = LASSO_all_genes(EUR_protein_genes, EUR_ge_regressed, "EUR")
    with open("./project_data/results/lasso_results.pkl", "wb") as f:
        pickle.dump(EUR_results, f)
    print("Results saved in ./project_data/results/lasso_results.pkl")

    YRI_results = LASSO_all_genes(YRI_protein_genes, YRI_ge_regressed, "YRI")
    with open("./project_data/results/lasso_results.pkl", "wb") as f:
        pickle.dump(YRI_results, f)
    print("Results saved in ./project_data/results/lasso_results.pkl")

# %%
