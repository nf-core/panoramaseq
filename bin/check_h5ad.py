#!/usr/bin/env python3
"""
check_h5ad_structure.py

Load an .h5ad file and verify its contents:
- Dimensions of X match obs and var
- obs contains 'batch'
- var index has no duplicates
- obsm contains 'spatial' with shape (n_obs,2) and checks for missing
"""

import argparse
import sys
import anndata as ad
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Check structure of an h5ad file")
    parser.add_argument("input", help="Path to the .h5ad file to check")
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    adata = ad.read_h5ad(args.input)

    # Check dimensions
    n_obs, n_vars = adata.n_obs, adata.n_vars
    print(f"adata.X shape: {adata.X.shape} (expected {n_obs} x {n_vars})")
    validation_errors = []
    
    if adata.X.shape != (n_obs, n_vars):
        error_msg = "ERROR: X shape mismatch!"
        print(error_msg)
        validation_errors.append(error_msg)

    # Check obs and var
    print(f"obs columns: {list(adata.obs.columns)}")
    if 'batch' not in adata.obs and 'sample' not in adata.obs:
        error_msg = "ERROR: Neither 'batch' nor 'sample' column found in obs."
        print(error_msg)
        validation_errors.append(error_msg)
    
    print(f"var index name: {adata.var.index.name}, number of genes: {n_vars}")
    dup_genes = adata.var.index.duplicated().sum()
    print(f"Duplicate genes in var index: {dup_genes}")
    if dup_genes > 0:
        error_msg = f"ERROR: Found {dup_genes} duplicate genes in var index."
        print(error_msg)
        validation_errors.append(error_msg)

    # Check spatial coords
    print(f"obsm keys: {list(adata.obsm.keys())}")
    if 'spatial' in adata.obsm:
        spatial = adata.obsm['spatial']
        print(f"Spatial coords shape: {spatial.shape} (expected {n_obs} x 2)")
        if spatial.shape != (n_obs, 2):
            error_msg = "ERROR: spatial shape mismatch!"
            print(error_msg)
            validation_errors.append(error_msg)
        n_missing = np.isnan(spatial).sum()
        print(f"Missing values in spatial coords: {n_missing}")
        if n_missing > n_obs * 0.5:  # More than 50% missing coordinates
            error_msg = f"ERROR: Too many missing spatial coordinates ({n_missing}/{n_obs*2})"
            print(error_msg)
            validation_errors.append(error_msg)
    else:
        error_msg = "ERROR: 'spatial' not found in obsm."
        print(error_msg)
        validation_errors.append(error_msg)

    # Final validation result
    if validation_errors:
        print(f"\n❌ VALIDATION FAILED with {len(validation_errors)} errors:")
        for error in validation_errors:
            print(f"  - {error}")
        print("H5AD file is invalid and will be removed.")
        sys.exit(1)
    else:
        print("\n✅ VALIDATION PASSED: H5AD file structure is correct.")
        print("All checks completed successfully.")


if __name__ == '__main__':
    main()
