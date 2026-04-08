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

    # Check if dataframe is empty or has no data rows
    if df.empty or len(df) == 0:
        raise ValueError(
            f"Input file '{path}' is empty or contains no data rows. "
            f"This typically happens when the GTF file has no overlapping features with aligned reads. "
            f"Please check that your GTF file contains relevant features for your data."
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
    empty_files = []
    for fp in args.inputs:
        sample_name = os.path.splitext(os.path.basename(fp))[0]
        try:
            adata_list.append(load_tsv_gz(fp, sample_name, coords_df))
        except ValueError as e:
            empty_files.append((fp, str(e)))
            print(f"WARNING: Skipping {fp}: {e}")

    # If all files are empty, exit with error
    if len(adata_list) == 0:
        print("\nERROR: All input files are empty or contain no data.")
        print("\nEmpty files encountered:")
        for fp, err in empty_files:
            print(f"  - {fp}: {err}")
        print("\nPossible causes:")
        print("  1. GTF file contains no features that overlap with aligned reads")
        print("  2. Test data is too small (e.g., using test profile with truncated GTF)")
        print("  3. Barcode calling filtered out all reads")
        print("  4. No UMIs were counted in the previous step")
        print("\nSuggestions:")
        print("  - Use a complete GTF file with full gene annotations")
        print("  - Check UMICOUNT logs for warnings")
        print("  - Verify that your barcode calling step produced valid output")
        raise SystemExit(1)

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
