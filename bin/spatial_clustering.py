#!/usr/bin/env python3
"""
Perform unsupervised clustering on a Panorama-seq h5ad and plot a spatial
heatmap with spots coloured by Leiden cluster.

Optional spatial binning (--bin-size N) aggregates nearby barcodes into
N×N pixel grid cells before clustering, reducing noise and memory use.
This mirrors the Visium HD bin2/bin8/bin16 approach.

Standard scanpy workflow:
  load → [bin] → filter → normalise → log1p → HVG → PCA → neighbours → Leiden → spatial plot

Usage:
    python3 spatial_clustering.py \
        --h5ad   <results>/h5ad/Sample.h5ad \
        --coords <results>/barcode_remapping/Sample_coords_synthetic.csv \
        --output cluster_map.png \
        [--bin-size 50] \
        [--resolution 0.5] \
        [--n-hvgs 2000] \
        [--min-counts 10] \
        [--n-pcs 30]
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_adata(h5ad_path: Path, coords_path: Path) -> ad.AnnData:
    print("Loading h5ad …")
    adata = ad.read_h5ad(h5ad_path)

    # Use gene symbols as var_names when index is Ensembl IDs
    if adata.var_names[0].startswith("ENSMUSG"):
        adata.var_names = adata.var["gene_name"].values.astype(str)
        adata.var_names_make_unique()

    # Attach spatial coordinates from the CSV
    print("Loading coordinates …")
    coords = pd.read_csv(coords_path, index_col=0)
    coords.index = coords.index.astype(str)
    adata.obs.index = adata.obs.index.astype(str)

    common = adata.obs.index.intersection(coords.index)
    if len(common) == 0:
        sys.exit("ERROR: No overlapping barcodes between h5ad and coordinate file.")

    adata = adata[common].copy()
    adata.obs["x"] = coords.loc[common, "x"].values.astype(float)
    adata.obs["y"] = coords.loc[common, "y"].values.astype(float)

    print(f"  {adata.n_obs:,} barcodes × {adata.n_vars:,} genes")
    n_pos = len(adata.obs[["x", "y"]].drop_duplicates())
    print(f"  Unique spatial positions: {n_pos:,}")
    return adata


def bin_adata(adata: ad.AnnData, bin_size: int) -> ad.AnnData:
    """
    Aggregate barcodes into spatial grid bins of size bin_size × bin_size pixels.
    Each bin becomes one observation; UMI counts are summed across barcodes in
    the same bin. The bin's (x, y) is the centre of the grid cell.
    """
    print(f"Binning into {bin_size}×{bin_size} pixel grid cells …")
    x = adata.obs["x"].values
    y = adata.obs["y"].values

    # Assign each barcode to a grid cell
    bin_x = (x // bin_size).astype(int)
    bin_y = (y // bin_size).astype(int)
    bin_id = [f"{bx}_{by}" for bx, by in zip(bin_x, bin_y)]
    adata.obs["bin_id"] = bin_id

    unique_bins = list(dict.fromkeys(bin_id))   # preserve order, deduplicate
    bin_index = {b: i for i, b in enumerate(unique_bins)}
    n_bins = len(unique_bins)

    # Build a barcode→bin mapping matrix and aggregate counts
    row_idx = np.array([bin_index[b] for b in bin_id], dtype=np.int32)
    col_idx = np.arange(adata.n_obs, dtype=np.int32)
    data    = np.ones(adata.n_obs, dtype=np.float32)
    agg_mat = csr_matrix((data, (row_idx, col_idx)),
                         shape=(n_bins, adata.n_obs))

    X_orig = adata.X if not issparse(adata.X) else adata.X
    X_binned = agg_mat @ X_orig          # (n_bins × n_genes), sparse or dense

    # Compute bin centre coordinates
    centres = (pd.DataFrame({"bin_id": bin_id, "x": x, "y": y})
               .groupby("bin_id")[["x", "y"]].mean())

    obs_df = pd.DataFrame({
        "x": centres.loc[unique_bins, "x"].values,
        "y": centres.loc[unique_bins, "y"].values,
    }, index=unique_bins)

    binned = ad.AnnData(
        X=csr_matrix(X_binned),
        obs=obs_df,
        var=adata.var.copy(),
    )

    n_counts = np.asarray(binned.X.sum(axis=1)).flatten()
    binned.obs["n_counts"] = n_counts

    print(f"  {adata.n_obs:,} barcodes → {n_bins:,} bins  "
          f"(mean {adata.n_obs / n_bins:.1f} barcodes/bin)")
    return binned


def run_clustering(adata: ad.AnnData,
                   min_counts: int,
                   n_hvgs: int,
                   n_pcs: int,
                   resolution: float) -> ad.AnnData:

    # --- Filter ---
    print(f"Filtering barcodes with < {min_counts} counts …")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After filter: {adata.n_obs:,} barcodes × {adata.n_vars:,} genes")

    # --- Normalise & log ---
    print("Normalising …")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --- Highly variable genes ---
    n_hvgs = min(n_hvgs, adata.n_vars)
    print(f"Selecting {n_hvgs} highly variable genes …")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()

    # --- PCA ---
    n_pcs = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    print(f"PCA ({n_pcs} PCs) …")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)

    # --- Neighbours & UMAP ---
    print("Computing neighbours …")
    sc.pp.neighbors(adata, n_pcs=n_pcs)

    # --- Leiden clustering ---
    print(f"Leiden clustering (resolution={resolution}) …")
    sc.tl.leiden(adata, resolution=resolution)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  {n_clusters} clusters found")

    return adata


def plot_spatial_clusters(adata: ad.AnnData, output: Path,
                          spot_size: float | None = None,
                          bin_size: int | None = None) -> None:
    clusters = adata.obs["leiden"].astype(str)
    unique_clusters = sorted(clusters.unique(), key=lambda x: int(x))
    n_clusters = len(unique_clusters)

    # Colour palette — tab20 for up to 20, then extend with tab20b/c
    cmap20  = plt.colormaps["tab20"].resampled(20)
    cmap20b = plt.colormaps["tab20b"].resampled(20)
    cmap20c = plt.colormaps["tab20c"].resampled(20)
    all_colors = ([cmap20(i)  for i in range(20)] +
                  [cmap20b(i) for i in range(20)] +
                  [cmap20c(i) for i in range(20)])
    color_map = {c: all_colors[i % len(all_colors)]
                 for i, c in enumerate(unique_clusters)}

    x = adata.obs["x"].values
    y = adata.obs["y"].values

    # Auto-size markers to fill the tissue area
    if spot_size is None:
        x_range = x.max() - x.min() if np.ptp(x) > 0 else 1
        y_range = y.max() - y.min() if np.ptp(y) > 0 else 1
        area_per_pt = (7 * 7 * 72**2) / max(adata.n_obs, 1)
        spot_size = min(area_per_pt * 0.5, 300)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                             gridspec_kw={"width_ratios": [3, 1]})
    ax, ax_leg = axes

    # Plot each cluster
    for cl in unique_clusters:
        mask = clusters == cl
        ax.scatter(x[mask], y[mask],
                   c=[color_map[cl]],
                   s=spot_size,
                   linewidths=0,
                   rasterized=True,
                   label=f"Cluster {cl}  (n={mask.sum():,})")

    ax.set_aspect("equal")
    ax.invert_yaxis()
    bin_label = f"  [bin {bin_size}px]" if bin_size else ""
    ax.set_title(f"Spatial clusters (Leiden){bin_label}", fontsize=13, fontweight="bold")
    ax.set_xlabel("x coordinate")
    ax.set_ylabel("y coordinate")
    ax.tick_params(labelsize=8)

    # Legend in second panel
    ax_leg.axis("off")
    handles = [mpatches.Patch(color=color_map[c], label=f"Cluster {c}  (n={( clusters == c).sum():,})")
               for c in unique_clusters]
    ax_leg.legend(handles=handles, loc="center left",
                  fontsize=8, frameon=False,
                  title=f"{n_clusters} clusters", title_fontsize=9)

    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad",       required=True, type=Path)
    parser.add_argument("--coords",     required=True, type=Path,
                        help="Barcode coordinate CSV (barcode, x, y)")
    parser.add_argument("--output",     default=Path("spatial_clusters.png"), type=Path)
    parser.add_argument("--bin-size",   default=None, type=int,
                        help="Aggregate barcodes into N×N pixel bins before clustering "
                             "(e.g. 50). Omit to cluster individual barcodes.")
    parser.add_argument("--resolution", default=0.5,  type=float,
                        help="Leiden resolution (default: 0.5)")
    parser.add_argument("--n-hvgs",     default=2000, type=int,
                        help="Number of highly variable genes (default: 2000)")
    parser.add_argument("--min-counts", default=10,   type=int,
                        help="Min UMI counts per barcode to keep (default: 10)")
    parser.add_argument("--n-pcs",      default=30,   type=int,
                        help="Number of PCA components (default: 30)")
    parser.add_argument("--spot-size",  default=None, type=float,
                        help="Scatter marker size in points² (auto if omitted)")
    parser.add_argument("--save-h5ad",  default=None, type=Path,
                        help="Optionally save the processed AnnData with cluster labels")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    adata = load_adata(args.h5ad, args.coords)

    if args.bin_size:
        adata = bin_adata(adata, args.bin_size)

    adata = run_clustering(adata,
                           min_counts=args.min_counts,
                           n_hvgs=args.n_hvgs,
                           n_pcs=args.n_pcs,
                           resolution=args.resolution)

    plot_spatial_clusters(adata, args.output,
                          spot_size=args.spot_size,
                          bin_size=args.bin_size)

    if args.save_h5ad:
        adata.write_h5ad(args.save_h5ad)
        print(f"Processed AnnData saved: {args.save_h5ad}")


if __name__ == "__main__":
    main()
