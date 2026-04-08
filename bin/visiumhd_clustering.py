#!/usr/bin/env python3
"""
Spatial clustering from a Visium HD feature_slice.h5 (Space Ranger output).

Reads count data directly from the feature_slice.h5, assembles a bins×genes
AnnData at the requested bin resolution (8 µm or 16 µm), runs the standard
scanpy clustering workflow, and plots a spatial heatmap coloured by cluster.

The script can also plot the pre-computed Space Ranger clusters stored inside
the h5 file (--use-spaceranger-clusters).

Usage:
    python3 visiumhd_clustering.py \
        --h5    Visium_HD_Mouse_Brain_Fresh_Frozen_feature_slice.h5 \
        --bin   016um \
        --output clusters_016um.png \
        [--resolution 0.5] \
        [--n-hvgs 3000] \
        [--min-counts 50] \
        [--n-pcs 30] \
        [--use-spaceranger-clusters]
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import coo_matrix, csr_matrix

# µm per pixel in the Visium HD 2µm grid
UM_PER_PIX = 2.0

VALID_BINS = ("008um", "016um")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_feature_slice(h5_path: Path, bin_res: str) -> ad.AnnData:
    """
    Parse the feature_slice.h5 and return an AnnData where:
      - obs  = spatial bins  (index = "row_col", columns = x_um, y_um)
      - var  = genes         (index = gene_name, columns = gene_id)
      - X    = UMI count matrix (bins × genes, sparse)

    The mask stores a binary sparse matrix: each nonzero (row, col) entry is
    a bin that is in tissue. row/col are coordinates in the *bin grid*, not
    in the 2µm pixel grid. We derive pixel→bin mapping by integer division:
      bin_row = pix_row // pix_per_bin
      bin_col = pix_col // pix_per_bin
    """
    if bin_res not in VALID_BINS:
        sys.exit(f"ERROR: --bin must be one of {VALID_BINS}, got '{bin_res}'")

    mask_key = f"square_{bin_res}"

    # pixels per bin edge (pixel grid is 2 µm/pixel)
    bin_um = int(bin_res.replace("um", ""))
    pix_per_bin = bin_um // int(UM_PER_PIX)   # e.g. 016um → 8 px, 008um → 4 px

    with h5py.File(h5_path, "r") as f:
        # ---- Gene metadata ----
        gene_names = np.array([g.decode() for g in f["features/name"][:]])
        gene_ids   = np.array([g.decode() for g in f["features/id"][:]])
        n_genes    = len(gene_names)

        # ---- Bin mask: (bin_row, bin_col) pairs that are in tissue ----
        if mask_key not in f["masks"]:
            sys.exit(f"ERROR: mask '{mask_key}' not found in h5. "
                     f"Available: {list(f['masks'].keys())}")

        mask      = f[f"masks/{mask_key}"]
        bin_rows  = mask["row"][:].astype(np.int32)   # bin-grid row
        bin_cols  = mask["col"][:].astype(np.int32)   # bin-grid col
        # data is all 1s (binary mask) — not used as bin IDs

        # Deduplicate and assign a sequential bin index
        bin_rc     = list(zip(bin_rows.tolist(), bin_cols.tolist()))
        unique_rc  = list(dict.fromkeys(bin_rc))          # ordered dedup
        bin_to_idx = {rc: i for i, rc in enumerate(unique_rc)}
        n_bins     = len(unique_rc)

        print(f"Bin resolution : {bin_res}  ({pix_per_bin} px/bin)")
        print(f"Tissue bins    : {n_bins:,}")

        # Bin coordinates in µm (centre of each bin cell)
        bin_y_um = np.array([(r + 0.5) * bin_um for r, c in unique_rc],
                            dtype=np.float32)
        bin_x_um = np.array([(c + 0.5) * bin_um for r, c in unique_rc],
                            dtype=np.float32)

        # Fast 2D lookup: bin_grid[bin_row, bin_col] → bin_index (or -1)
        max_br = int(bin_rows.max()) + 1
        max_bc = int(bin_cols.max()) + 1
        bin_grid = np.full((max_br, max_bc), -1, dtype=np.int32)
        for idx, (r, c) in enumerate(unique_rc):
            bin_grid[r, c] = idx

        # ---- Read feature slices and accumulate counts ----
        fs_keys        = list(f["feature_slices"].keys())
        valid_gene_idx = sorted([int(k) for k in fs_keys])

        print(f"Reading {len(valid_gene_idx):,} / {n_genes:,} expressed genes …")

        rows_out, cols_out, data_out = [], [], []

        for i, gene_idx in enumerate(valid_gene_idx):
            if i % 2000 == 0:
                print(f"  gene {i:,} / {len(valid_gene_idx):,} …")

            fs = f[f"feature_slices/{gene_idx}"]
            pix_r  = fs["row"][:].astype(np.int32)
            pix_c  = fs["col"][:].astype(np.int32)
            counts = fs["data"][:].astype(np.float32)

            # Pixel → bin-grid coords (vectorised)
            br = pix_r // pix_per_bin
            bc = pix_c // pix_per_bin

            # Clip and look up in bin_grid
            in_bounds = (br < max_br) & (bc < max_bc)
            br = br[in_bounds]; bc = bc[in_bounds]; counts = counts[in_bounds]

            bidx = bin_grid[br, bc]
            keep = bidx >= 0
            if not keep.any():
                continue
            br = br[keep]; bc = bc[keep]
            bidx = bidx[keep]; counts = counts[keep]

            # Aggregate multiple pixels in the same bin
            gene_vec = np.zeros(n_bins, dtype=np.float32)
            np.add.at(gene_vec, bidx, counts)

            nz = gene_vec.nonzero()[0]
            rows_out.append(nz)
            cols_out.append(np.full(len(nz), gene_idx, dtype=np.int32))
            data_out.append(gene_vec[nz])

        print("Assembling count matrix …")
        all_rows = np.array(rows_out, dtype=np.int32)
        all_cols = np.array(cols_out, dtype=np.int32)
        all_data = np.array(data_out, dtype=np.float32)

        X = coo_matrix(
            (all_data, (all_rows, all_cols)),
            shape=(n_bins, n_genes),
        ).tocsr()

    # ---- Build AnnData ----
    bin_labels = [f"{r}_{c}" for r, c in unique_rc]

    obs_df = pd.DataFrame({
        "x_um": bin_x_um,
        "y_um": bin_y_um,
    }, index=bin_labels)

    var_df = pd.DataFrame({"gene_id": gene_ids}, index=gene_names)
    var_df.index.name = None

    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    adata.var_names_make_unique()

    n_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    adata.obs["n_counts"] = n_counts
    adata = adata[adata.obs["n_counts"] > 0].copy()

    print(f"AnnData: {adata.n_obs:,} bins × {adata.n_vars:,} genes")
    return adata


def load_spaceranger_clusters(h5_path: Path, bin_res: str,
                              n_bins: int) -> pd.Series | None:
    """
    Load graph-based clusters pre-computed by Space Ranger.
    Returns a Series indexed by bin_id (str), or None if not found.
    """
    key = f"square_{bin_res}_gene_expression_graphclust"
    with h5py.File(h5_path, "r") as f:
        clust_keys = list(f["secondary_analysis/clustering"].keys())
        if key not in clust_keys:
            print(f"WARNING: Space Ranger cluster key '{key}' not found. "
                  f"Available: {clust_keys}")
            return None

        cl = f[f"secondary_analysis/clustering/{key}"]
        rows = cl["row"][:]    # bin_id
        data = cl["data"][:]   # cluster label
    sr = pd.Series(data.astype(str), index=rows.astype(str),
                   name="spaceranger_cluster")
    return sr


# ---------------------------------------------------------------------------
# Clustering (reuse same scanpy workflow)
# ---------------------------------------------------------------------------

def run_clustering(adata: ad.AnnData,
                   min_counts: int,
                   n_hvgs: int,
                   n_pcs: int,
                   resolution: float) -> ad.AnnData:

    print(f"Filtering bins with < {min_counts} counts …")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After filter: {adata.n_obs:,} bins × {adata.n_vars:,} genes")

    print("Normalising …")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_hvgs = min(n_hvgs, adata.n_vars)
    print(f"Selecting {n_hvgs} HVGs …")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()

    n_pcs = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    print(f"PCA ({n_pcs} PCs) …")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)

    print("Computing neighbours …")
    sc.pp.neighbors(adata, n_pcs=n_pcs)

    print(f"Leiden clustering (resolution={resolution}) …")
    sc.tl.leiden(adata, resolution=resolution,
                 flavor="igraph", n_iterations=2, directed=False)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  {n_clusters} clusters found")
    return adata


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _make_color_map(unique_labels):
    cmap20  = plt.colormaps["tab20"].resampled(20)
    cmap20b = plt.colormaps["tab20b"].resampled(20)
    cmap20c = plt.colormaps["tab20c"].resampled(20)
    colors  = ([cmap20(i)  for i in range(20)] +
               [cmap20b(i) for i in range(20)] +
               [cmap20c(i) for i in range(20)])
    return {c: colors[i % len(colors)] for i, c in enumerate(unique_labels)}


def plot_clusters(adata: ad.AnnData, cluster_col: str,
                  output: Path, title: str,
                  spot_size: float | None = None) -> None:

    labels = adata.obs[cluster_col].astype(str)
    try:
        unique = sorted(labels.unique(), key=lambda x: int(x))
    except ValueError:
        unique = sorted(labels.unique())

    color_map = _make_color_map(unique)

    x = adata.obs["x_um"].values
    y = adata.obs["y_um"].values

    if spot_size is None:
        area_per_pt = (8 * 8 * 72**2) / max(adata.n_obs, 1)
        spot_size = min(area_per_pt * 0.4, 200)

    fig, (ax, ax_leg) = plt.subplots(
        1, 2, figsize=(15, 8),
        gridspec_kw={"width_ratios": [3, 1]})

    for cl in unique:
        mask = labels == cl
        ax.scatter(x[mask], y[mask],
                   c=[color_map[cl]], s=spot_size,
                   linewidths=0, rasterized=True)

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.tick_params(labelsize=8)

    ax_leg.axis("off")
    handles = [mpatches.Patch(color=color_map[c],
                              label=f"Cluster {c}  (n={(labels == c).sum():,})")
               for c in unique]
    ax_leg.legend(handles=handles, loc="center left", fontsize=7,
                  frameon=False,
                  title=f"{len(unique)} clusters", title_fontsize=9)

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
    parser.add_argument("--h5",     required=True, type=Path,
                        help="Path to feature_slice.h5")
    parser.add_argument("--bin",    default="016um", choices=VALID_BINS,
                        dest="bin_res",
                        help="Bin resolution to use (default: 016um)")
    parser.add_argument("--output", default=Path("visiumhd_clusters.png"), type=Path)
    parser.add_argument("--resolution", default=0.5,  type=float)
    parser.add_argument("--n-hvgs",     default=3000, type=int)
    parser.add_argument("--min-counts", default=50,   type=int)
    parser.add_argument("--n-pcs",      default=30,   type=int)
    parser.add_argument("--spot-size",  default=None, type=float)
    parser.add_argument("--save-h5ad",  default=None, type=Path)
    parser.add_argument("--use-spaceranger-clusters", action="store_true",
                        help="Also plot the Space Ranger pre-computed clusters")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # ---- Load data ----
    adata = read_feature_slice(args.h5, args.bin_res)

    # ---- Cluster ----
    adata = run_clustering(adata,
                           min_counts=args.min_counts,
                           n_hvgs=args.n_hvgs,
                           n_pcs=args.n_pcs,
                           resolution=args.resolution)

    # ---- Plot Leiden clusters ----
    title = (f"Leiden clusters  [{args.bin_res} bins, "
             f"res={args.resolution}]")
    plot_clusters(adata, "leiden", args.output, title,
                  spot_size=args.spot_size)

    # ---- Optionally plot Space Ranger clusters ----
    if args.use_spaceranger_clusters:
        sr_clusters = load_spaceranger_clusters(
            args.h5, args.bin_res, adata.n_obs)
        if sr_clusters is not None:
            # Align to filtered adata
            common = adata.obs.index.intersection(sr_clusters.index)
            adata_sr = adata[common].copy()
            adata_sr.obs["sr_cluster"] = sr_clusters.loc[common].values
            sr_out = args.output.with_stem(args.output.stem + "_spaceranger")
            plot_clusters(adata_sr, "sr_cluster", sr_out,
                          f"Space Ranger clusters  [{args.bin_res}]",
                          spot_size=args.spot_size)

    # ---- Save processed AnnData ----
    if args.save_h5ad:
        adata.write_h5ad(args.save_h5ad)
        print(f"Saved AnnData: {args.save_h5ad}")


if __name__ == "__main__":
    main()
