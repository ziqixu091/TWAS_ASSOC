import subprocess
import pandas as pd
import numpy as np
import random
from helper_fns import *

from concurrent.futures import ProcessPoolExecutor
import concurrent

from tqdm import tqdm

random.seed(42)

# Load gene definitions
EUR_ge_regressed, YRI_ge_regressed, EUR_protein_genes, YRI_protein_genes = load_data()

os.chdir("/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific")

# Function to extract SNPs for a gene and split into train/test
# def process_gene(gene, chr, start, end, ancestry="EUR", plink_path='plink'):
def process_gene(row, ancestry="EUR", plink_path='plink'):
    gene = row['gene_id']
    chr = row['chr']
    start = row['start']
    end = row['end']

    if start < 0 or start > end:
        print(f"Invalid start/end for {gene}")
        return
    
    os.makedirs(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene}', exist_ok=True)
    os.chdir(f'/new-stg/home/banghua/TWAS_ASSOC/project_data/geno_gene_specific/{ancestry}/{gene}')


    base_name = f'{gene}_chr{chr}_{start}_{end}'
    if ancestry == "EUR":
        bfile_path = "/new-stg/home/banghua/TWAS_ASSOC/project_data/geno/EUR/GEUVADIS_EUR_chr" + str(chr)
    else:
        bfile_path = "/new-stg/home/banghua/TWAS_ASSOC/project_data/geno/YRI/GEUVADIS_YRI_chr" + str(chr)

    if ancestry == "EUR":
        gene_expression = EUR_ge_regressed
    else:
        gene_expression = YRI_ge_regressed

    gene_exp_one_gene = gene_expression[[gene]]

    # Step 2: Extract gene-specific SNPs
    try:
        subprocess.run([plink_path, '--bfile', bfile_path, '--chr', str(chr), 
                        '--from-bp', str(start), '--to-bp', str(end), '--allow-no-sex',
                        '--make-bed', '--out', base_name])
    

        # Step 3: Add gene expression data
        fam = pd.read_csv(f"./{base_name}.fam", sep='\s+', header=None)
        all_individuals = fam[[1]].values.flatten().tolist()
        fam.index = all_individuals
        gene_exp_one_gene = gene_exp_one_gene.loc[all_individuals]
        # Add gene expression to fam
        fam[[5]] = gene_exp_one_gene
        fam.to_csv(f"./{base_name}.fam", sep=' ', index=False, header=False)

        # Count the number of individuals
        with open(f'{base_name}.fam') as f:
            num_individuals = sum(1 for line in f)
        
        # Calculate split
        train_num = int(num_individuals * 0.8)
        
        # Shuffle individuals
        individuals = pd.read_csv(f'{base_name}.fam', sep='\s+', header=None)
        shuffled = individuals.sample(frac=1).reset_index(drop=True)
        
        # Split into train/test
        train = shuffled.head(train_num)
        test = shuffled.tail(num_individuals - train_num)
        
        train[[0, 1]].to_csv(f'{base_name}_train.txt', sep=' ', index=False, header=False)
        test[[0, 1]].to_csv(f'{base_name}_test.txt', sep=' ', index=False, header=False)
        
        # Generate train/test datasets
        subprocess.run([plink_path, '--bfile', base_name, '--keep', f'{base_name}_train.txt', '--allow-no-sex',
                        '--make-bed', '--out', f'{base_name}_train'])
        subprocess.run([plink_path, '--bfile', base_name, '--keep', f'{base_name}_test.txt', '--allow-no-sex',
                        '--make-bed', '--out', f'{base_name}_test'])
    except:
        print(f"Error processing {gene}")
        # Remove the directory (with files) if there was an error
        os.chdir('..')
        subprocess.run(['rm', '-r', gene])


if __name__ == "__main__":
    # Use ProcessPoolExecutor to execute the tasks in parallel
    with ProcessPoolExecutor() as executor:
        # Using list comprehension to create and immediately start all tasks
        # tqdm is used to display progress
        futures = [executor.submit(process_gene, row, "EUR") for index, row in EUR_protein_genes.iterrows()]
        # Ensure all tasks are complete before moving on
        results = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            results.append(future.result())

        
    with ProcessPoolExecutor() as executor:
        # Using list comprehension to create and immediately start all tasks
        # tqdm is used to display progress
        futures = [executor.submit(process_gene, row, "YRI") for index, row in YRI_protein_genes.iterrows()]
        # Ensure all tasks are complete before moving on
        results = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            results.append(future.result())