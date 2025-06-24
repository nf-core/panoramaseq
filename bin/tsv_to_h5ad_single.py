#!/usr/bin/env python3
"""
tsv_to_h5ad_single.py

Convert a single gzipped 3-column count TSV (gene, cell, count)
into an AnnData (.h5ad), and attach spatial coordinates from
a separate CSV (cell, x, y).
"""

import os
import argparse
import pandas as pd
from scipy import sparse
import anndata as ad

def load_single_tsv_gz(path, sample_name, coords_df):
    """
    Load a 3-column TSV (gene, cell, count) into an AnnData,
    and attach spatial coords from coords_df.
    """
    # 1) Read counts
    df = pd.read_csv(
        path, sep="\t", compression="gzip",
        usecols=["gene", "cell", "count"]
    )

    # 2) Build gene & cell indices
    genes = pd.Index(df["gene"].unique(), name="gene")
    cells = pd.Index(df["cell"].unique(), name="cell")

    # 3) Map to integer arrays
    gene_idx = genes.get_indexer(df["gene"])
    cell_idx = cells.get_indexer(df["cell"])

    # 4) Build sparse count matrix: rows=cells, cols=genes
    X = sparse.csr_matrix(
        (df["count"].values, (cell_idx, gene_idx)),
        shape=(len(cells), len(genes))
    )

    # 5) Make AnnData with obs (cells) and var (genes)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes)
    )

    # 6) Attach spatial coords (will align on the original cell names)
    #    any missing coords will become NaN
    sub = coords_df.reindex(cells)  # index=cell
    # store in adata.obsm["spatial"] as an (n_obs × 2) array
    adata.obsm["spatial"] = sub[["x", "y"]].to_numpy()

    # 7) Annotate sample
    adata.obs["sample"] = sample_name

    return adata

def main():
    p = argparse.ArgumentParser(
        description="Convert single gzipped count TSV into .h5ad, with spatial coords"
    )
    p.add_argument(
        "input", 
        help="Input gzipped count TSV (*.tsv.gz) with columns gene, cell, count"
    )
    p.add_argument(
        "-c", "--coords", required=True,
        help="CSV file (cell, x, y) mapping each cell barcode to spatial coordinates"
    )
    p.add_argument(
        "-o", "--output", required=True,
        help="Output HDF5 AnnData file (e.g. sample.h5ad)"
    )
    p.add_argument(
        "-s", "--sample-name",
        help="Sample name to use in metadata (defaults to input filename without extension)"
    )
    args = p.parse_args()

    # Load coords once
    coords_df = pd.read_csv(args.coords, usecols=["cell", "x", "y"])
    coords_df.set_index("cell", inplace=True)

    # Determine sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = os.path.splitext(os.path.basename(args.input))[0]

    # Process the sample
    adata = load_single_tsv_gz(args.input, sample_name, coords_df)

    # Write to .h5ad
    adata.write_h5ad(args.output)
    print(
        f"Written AnnData with {adata.n_obs} cells, "
        f"{adata.n_vars} genes → {args.output}"
    )

if __name__ == "__main__":
    main()
