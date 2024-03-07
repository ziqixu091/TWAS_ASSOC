import subprocess
import pandas as pd
import numpy as np
import random
from helper_fns import *
random.seed(42)

# Function to extract SNPs for a gene and split into train/test
def process_gene(gene, chr, start, end, full_y_df, test_size=0.2, ancestry='EUR', plink_path='/new-stg/home/banghua/anaconda3/envs/CSE_284/bin/plink'):
    # base_name = f'{gene}_chr{chr}_{start}_{end}'
    if not os.path.exists(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene}'):
        os.makedirs(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene}')
    
    b_file_path = "/new-stg/home/banghua/TWAS_ASSOC/project_data/geno/" + ancestry + "/GEUVADIS_EUR_chr" + str(chr)
    out_b_file_path = f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene}/{gene}'
    # Step 2: Extract gene-specific SNPs
    subprocess.run([plink_path, '--bfile', b_file_path, '--chr', str(chr), 
                    '--from-bp', str(start), '--to-bp', str(end), 
                    '--make-bed', '--out', out_b_file_path])
    
    # Count the number of individuals
    with open(f'{out_b_file_path}.fam') as f:
        num_individuals = sum(1 for line in f)
    
    # Calculate split
    train_num = int(num_individuals * (1-float(test_size)))
    
    # Shuffle individuals
    individuals = pd.read_csv(f'{out_b_file_path}.fam', sep='\s+', header=None)
    shuffled = individuals.sample(frac=1).reset_index(drop=True)
    
    # Split into train/test
    train = shuffled.head(train_num)
    test = shuffled.tail(num_individuals - train_num)
    
    train[[0, 1]].to_csv(f'{out_b_file_path}_train.txt', sep='\t', index=False, header=False)
    test[[0, 1]].to_csv(f'{out_b_file_path}_test.txt', sep='\t', index=False, header=False)
    
    # Generate train/test datasets
    subprocess.run([plink_path, '--bfile', out_b_file_path, '--keep', f'{out_b_file_path}_train.txt', 
                    '--make-bed', '--out', f'{out_b_file_path}_train'])
    subprocess.run([plink_path, '--bfile', out_b_file_path, '--keep', f'{out_b_file_path}_test.txt', 
                    '--make-bed', '--out', f'{out_b_file_path}_test'])
    
    # Step 3: Generate gene-specific phenotype
    # Get gene expression
    gene_idx = full_y_df.columns.get_loc(gene)
    gene_expression = full_y_df.iloc[:, gene_idx]

    # Generate gene-specific phenotype
    train_phenotype = gene_expression[train[1]]
    test_phenotype = gene_expression[test[1]]

    # Save gene-specific phenotype
    train_phenotype.to_csv(f'{out_b_file_path}_train.pheno', sep='\t', header=False)
    test_phenotype.to_csv(f'{out_b_file_path}_test.pheno', sep='\t', header=False)


if __name__ == "__main__":
    # Load data
    EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes = load_data()

    for gene in EUR_protein_genes['gene_id']:
        if os.path.exists(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/EUR/{gene}'):
            continue
        gene_info = EUR_protein_genes[EUR_protein_genes['gene_id'] == gene]
        chr_num = gene_info['chr'].values[0]
        start = gene_info['start'].values[0]
        end = gene_info['end'].values[0]
        try:
            process_gene(gene, chr_num, start, end, EUR_ge_regressed, ancestry='EUR')
        except:
            print(f'Error processing gene {gene} in EUR.')
            os.rmdir(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/EUR/{gene}')
            continue

    print(f'Finished processing {len(EUR_protein_genes)} EUR genes.')

    for gene in YRI_protein_genes['gene_id']:
        if os.path.exists(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/YRI/{gene}'):
            continue
        gene_info = YRI_protein_genes[YRI_protein_genes['gene_id'] == gene]
        chr_num = gene_info['chr'].values[0]
        start = gene_info['start'].values[0]
        end = gene_info['end'].values[0]
        try:
            process_gene(gene, chr_num, start, end, YRI_ge_regressed, ancestry='YRI')
        except:
            print(f'Error processing gene {gene} in YRI.')
            os.rmdir(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/YRI/{gene}')
            continue

    print(f'Finished processing {len(YRI_protein_genes)} YRI genes.')
