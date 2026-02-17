#!/usr/bin/env python3

import argparse
import pandas as pd
import scipy.io as sio
import anndata
import numpy as np
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Convert STARsolo output to H5AD with spatial coordinates')
    parser.add_argument('--matrix', required=True, help='Path to matrix.mtx file')
    parser.add_argument('--barcodes', required=True, help='Path to barcodes.tsv file')
    parser.add_argument('--features', required=True, help='Path to features.tsv file')
    parser.add_argument('--coords', required=True, help='Path to barcode coordinates CSV file')
    parser.add_argument('--output', required=True, help='Output H5AD file path')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    print(f"=== Converting STARsolo Output to H5AD ===")
    print(f"Sample: {args.sample_id}")
    print(f"Matrix: {args.matrix}")
    print(f"Barcodes: {args.barcodes}")
    print(f"Features: {args.features}")
    print(f"Coordinates: {args.coords}")
    print(f"==========================================")
    
    # Read STARsolo output (Market Matrix format)
    print("Reading matrix...")
    matrix = sio.mmread(args.matrix).T.tocsr()
    
    print("Reading barcodes...")
    barcodes = pd.read_csv(args.barcodes, header=None, names=["barcode"], sep="\\t")
    
    print("Reading features...")
    features = pd.read_csv(args.features, sep="\\t", header=None)
    if features.shape[1] >= 2:
        features.columns = ["gene_id", "gene_name"] + [f"col{i}" for i in range(2, features.shape[1])]
        features = features[["gene_id", "gene_name"]]
    else:
        features.columns = ["gene_id"]
        features["gene_name"] = features["gene_id"]
    
    print(f"Matrix shape: {matrix.shape}")
    print(f"Number of barcodes: {len(barcodes)}")
    print(f"Number of features: {len(features)}")
    
    # Create AnnData object
    print("Creating AnnData object...")
    adata = anndata.AnnData(X=matrix, obs=barcodes, var=features)
    # Use gene_id as index (unique), keep gene_name as a column
    adata.var_names = adata.var["gene_id"]
    adata.obs_names = adata.obs["barcode"]
    
    # Add spatial coordinates from barcode file
    print("Adding spatial coordinates...")
    coords = pd.read_csv(args.coords)
    
    # Ensure barcode column exists
    if 'barcode' not in coords.columns:
        # Assume first column is barcode
        coords.columns = ['barcode'] + list(coords.columns[1:])
    
    # Set barcode as index in coords for merging
    coords = coords.set_index('barcode')
    
    # Merge coordinates with obs using index
    adata.obs = adata.obs.join(coords, how="left")
    
    # Check for spatial columns
    spatial_cols = [col for col in adata.obs.columns if col.lower() in ['x', 'y', 'x_coord', 'y_coord', 'row', 'col']]
    if spatial_cols:
        print(f"Found spatial columns: {spatial_cols}")
        adata.obs['has_spatial'] = ~adata.obs[spatial_cols].isnull().all(axis=1)
    else:
        print("Warning: No spatial coordinate columns found")
    
    # Add sample metadata
    adata.obs['sample_id'] = args.sample_id
    
    # Calculate basic QC metrics manually
    print("Calculating QC metrics...")
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata.var['total_counts'] = np.array(adata.X.sum(axis=0)).flatten()
    
    # Save to H5AD
    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    
    print("=== Conversion Complete ===")
    print(f"Total cells: {adata.n_obs}")
    print(f"Total genes: {adata.n_vars}")
    print(f"Total UMIs: {adata.X.sum():.0f}")
    print(f"Mean UMIs per cell: {adata.X.sum(axis=1).mean():.2f}")
    print(f"Median UMIs per cell: {pd.Series(adata.X.sum(axis=1).A1).median():.2f}")
    print("===========================")

if __name__ == "__main__":
    main()
