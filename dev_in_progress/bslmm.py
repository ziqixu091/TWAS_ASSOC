import subprocess
import pandas as pd
import numpy as np

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

    # Constructing the GEMMA command
    arg = f"{gemma_path} -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile {input} -bslmm {bv_type} -o {out}"

    # Execute the GEMMA command
    result = subprocess.run(arg, shell=True, capture_output=not sys_print)
    if not sys_print:
        print(result.stdout.decode())  # Optional: print GEMMA output for debugging.

    # Read the output parameter file
    try:
        eff = pd.read_table(f"{out}.param.txt", header=0, sep='\t')
    except FileNotFoundError:
        raise FileNotFoundError("GEMMA output file not found. Check GEMMA execution and output path.")

    # Initialize effect weights with NaN for all SNPs
    eff_wgt = pd.Series(np.nan, index=snp)

    # Match SNPs and assign weights
    for i, snp_id in enumerate(snp):
        if snp_id in eff['rs'].values:
            row = eff.loc[eff['rs'] == snp_id].iloc[0]
            eff_wgt.at[snp_id] = row['alpha'] + row['beta'] * row['gamma']

    return eff_wgt.values
