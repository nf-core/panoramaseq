#!/usr/bin/env python3
"""
tsv_to_h5ad.py

Convert one or more gzipped 3-column count TSVs (gene, cell, count)
into a merged AnnData (.h5ad), and attach spatial coordinates from
a separate CSV (cell, x, y). Each input will become one “batch”
with cell names prefixed by the filename.
"""

import os
import argparse
import pandas as pd
from scipy import sparse
import anndata as ad

def load_tsv_gz(path, prefix, coords_df):
    """
    Load a 3-column TSV (gene, cell, count) into an AnnData,
    using 'prefix' to tag the batch and make obs_names unique,
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

    # 7) Annotate batch and prefix obs_names
    adata.obs["batch"] = prefix
    adata.obs_names = [f"{prefix}_{bc}" for bc in adata.obs_names]

    return adata

def main():
    p = argparse.ArgumentParser(
        description="Merge gzipped count TSVs into one .h5ad, with spatial coords"
    )
    p.add_argument(
        "inputs", nargs="+",
        help="Input gzipped count TSVs (*.tsv.gz), each with columns gene, cell, count"
    )
    p.add_argument(
        "-c", "--coords", required=True,
        help="CSV file (cell, x, y) mapping each cell barcode to spatial coordinates"
    )
    p.add_argument(
        "-o", "--output", required=True,
        help="Output HDF5 AnnData file (e.g. merged_with_spatial.h5ad)"
    )
    args = p.parse_args()

    # Load coords once
    coords_df = pd.read_csv(args.coords, usecols=["cell", "x", "y"])
    coords_df.set_index("cell", inplace=True)

    # Process each sample
    adata_list = []
    for fp in args.inputs:
        sample_name = os.path.splitext(os.path.basename(fp))[0]
        adata_list.append(load_tsv_gz(fp, sample_name, coords_df))

    # Concatenate all samples (outer join on genes; obs already unique)
    combined = ad.concat(
        adata_list,
        join="outer",
        merge="unique",
        label="batch"
    )

    # Write to .h5ad
    combined.write_h5ad(args.output)
    print(
        f"Written merged AnnData with {combined.n_obs} cells, "
        f"{combined.n_vars} genes → {args.output}"
    )

if __name__ == "__main__":
    main()
