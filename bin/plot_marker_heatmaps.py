#!/usr/bin/env python3
"""
Plot per-spot spatial expression heatmaps for marker genes from Panorama-seq h5ad output.

Usage:
    python3 plot_marker_heatmaps.py \
        --h5ad   <results>/h5ad/Sample.h5ad \
        --coords <results>/barcode_remapping/Sample_coords_synthetic.csv \
        --outdir <output_directory>
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.sparse import issparse

# ---------------------------------------------------------------------------
# Marker gene definitions  (Ctip2 is an alias for Bcl11b — kept as note)
# ---------------------------------------------------------------------------
MARKERS = {
    "Cortex": ["Mbp", "Plp1"],
    "Hippocampus": ["Neurod6", "Zbtb20"],
    "Striatum": ["Foxp1", "Adora2a","Pdyn"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_adata(h5ad_path: Path, coords_path: Path) -> ad.AnnData:
    """Load h5ad and attach spatial coordinates from the CSV file."""
    adata = ad.read_h5ad(h5ad_path)

    # Use gene_name as var_names if the current index is Ensembl IDs
    if adata.var_names[0].startswith("ENSMUSG"):
        adata.var_names = adata.var["gene_name"].values.astype(str)
        adata.var_names_make_unique()

    # Load coordinates
    coords = pd.read_csv(coords_path, index_col=0)   # index = barcode
    coords.index = coords.index.astype(str)

    # Align to adata.obs
    adata.obs.index = adata.obs.index.astype(str)
    common = adata.obs.index.intersection(coords.index)
    if len(common) == 0:
        sys.exit("ERROR: No overlapping barcodes between h5ad and coordinate file.")

    adata = adata[common].copy()
    adata.obs["x"] = coords.loc[common, "x"].values
    adata.obs["y"] = coords.loc[common, "y"].values

    print(f"Loaded {adata.n_obs:,} barcodes x {adata.n_vars:,} genes")
    n_pos = len(adata.obs[["x", "y"]].drop_duplicates())
    print(f"Unique spatial positions: {n_pos}")
    return adata


def gene_expression_by_position(adata: ad.AnnData, gene: str) -> pd.DataFrame:
    """Return a DataFrame with columns [x, y, expr] aggregated (sum) per spot."""
    idx = adata.var_names.get_loc(gene)
    mat = adata.X[:, idx]
    if issparse(mat):
        expr = np.asarray(mat.todense()).flatten()
    else:
        expr = np.asarray(mat).flatten()

    df = pd.DataFrame({"x": adata.obs["x"].values,
                        "y": adata.obs["y"].values,
                        "expr": expr})
    agg = df.groupby(["x", "y"], as_index=False)["expr"].sum()
    return agg


def plot_gene(agg: pd.DataFrame, gene: str, region: str,
              outdir: Path, vmax_pct: float = 99.0) -> None:
    """Save a spatial scatter heatmap for a single gene."""
    fig, ax = plt.subplots(figsize=(6, 6))

    # Colour scale: cap at vmax_pct percentile to avoid outlier washout
    vmax = np.percentile(agg["expr"], vmax_pct) if agg["expr"].max() > 0 else 1.0
    vmax = max(vmax, 1.0)

    # Marker size: fill available area with circles so spots touch
    x_range = agg["x"].max() - agg["x"].min() if agg["x"].nunique() > 1 else 1
    y_range = agg["y"].max() - agg["y"].min() if agg["y"].nunique() > 1 else 1
    n_pts = len(agg)
    # Approximate: divide canvas area by number of points, convert to points²
    area_per_pt = (6 * 6 * 72**2) / max(n_pts, 1)   # figure 6×6 inch at 72 dpi
    marker_area = min(area_per_pt * 0.6, 400)

    # Plot zero-count spots as light grey background
    zero = agg[agg["expr"] == 0]
    if len(zero) > 0:
        ax.scatter(zero["x"], zero["y"], c="#d0d0d0", s=marker_area,
                   linewidths=0, rasterized=True, zorder=1)

    # Plot only spots with non-zero expression on top
    pos = agg[agg["expr"] > 0]
    sc = ax.scatter(
        pos["x"], pos["y"],
        c=pos["expr"],
        s=marker_area,
        cmap="Reds",
        vmin=0.5,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )

    if len(pos) > 0:
        cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("UMI count (sum per spot)", fontsize=9)
    else:
        ax.text(0.5, 0.5, "no counts detected", transform=ax.transAxes,
                ha="center", va="center", color="grey", fontsize=11)

    ax.set_aspect("equal")
    ax.invert_yaxis()          # image convention: y increases downward
    ax.set_title(f"{gene}  [{region}]", fontsize=13, fontweight="bold")
    ax.set_xlabel("x coordinate")
    ax.set_ylabel("y coordinate")
    ax.tick_params(labelsize=8)

    fig.tight_layout()
    out_file = outdir / f"{region}_{gene}.png"
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad",   required=True, type=Path, help="Path to .h5ad file")
    parser.add_argument("--coords", required=True, type=Path,
                        help="Path to barcode coordinate CSV (barcode, x, y)")
    parser.add_argument("--outdir", default=Path("marker_heatmaps"), type=Path,
                        help="Output directory for PNG files (default: marker_heatmaps/)")
    parser.add_argument("--vmax-pct", default=99.0, type=float,
                        help="Percentile for colour scale upper bound (default: 99)")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    adata = load_adata(args.h5ad, args.coords)

    # Check which markers are present
    all_markers = [g for genes in MARKERS.values() for g in genes]
    present = [g for g in all_markers if g in adata.var_names]
    absent  = [g for g in all_markers if g not in adata.var_names]
    if absent:
        print(f"WARNING: genes not found and will be skipped: {absent}")

    # Plot
    for region, genes in MARKERS.items():
        region_dir = args.outdir / region
        region_dir.mkdir(exist_ok=True)
        print(f"\n--- {region} ---")
        for gene in genes:
            if gene not in adata.var_names:
                print(f"  SKIP {gene} (not in dataset)")
                continue
            agg = gene_expression_by_position(adata, gene)
            total_expr = agg["expr"].sum()
            n_nonzero  = (agg["expr"] > 0).sum()
            print(f"  {gene}: total UMIs={total_expr:.0f}, non-zero spots={n_nonzero}/{len(agg)}")
            plot_gene(agg, gene, region, region_dir, vmax_pct=args.vmax_pct)

    print(f"\nDone. All figures saved to: {args.outdir}/")


if __name__ == "__main__":
    main()
