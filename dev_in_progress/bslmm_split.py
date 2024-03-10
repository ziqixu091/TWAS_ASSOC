import subprocess
import pandas as pd
import numpy as np
import random
from helper_fns import *

random.seed(42)

# Load gene definitions
gene_defs = pd.read_csv('gene_definitions.txt', sep='\t')

# Function to extract SNPs for a gene and split into train/test
def process_gene(gene, chr, start, end, plink_path='plink'):
    base_name = f'{gene}_chr{chr}_{start}_{end}'
    # Step 2: Extract gene-specific SNPs
    subprocess.run([plink_path, '--bfile', f'chr{chr}', '--chr', str(chr), 
                    '--from-bp', str(start), '--to-bp', str(end), 
                    '--make-bed', '--out', base_name])
    
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
    subprocess.run([plink_path, '--bfile', base_name, '--keep', f'{base_name}_train.txt', 
                    '--make-bed', '--out', f'{base_name}_train'])
    subprocess.run([plink_path, '--bfile', base_name, '--keep', f'{base_name}_test.txt', 
                    '--make-bed', '--out', f'{base_name}_test'])

# Iterate through each gene
for index, row in gene_defs.iterrows():
    process_gene(row['GeneName'], row['Chr'], row['Start'], row['End'])

