#!/usr/bin/env python3
"""
Plot spatial expression heatmaps for marker genes directly from a
Visium HD feature_slice.h5 file - no pipeline output required.

Usage:
    python3 plot_marker_heatmaps_visiumhd.py \
        --h5   <path>/feature_slice.h5 \
        --outdir real_visiumhd_heatmaps/
"""

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MARKERS = {
    "Cortex":      ["Reln", "Cux1", "Cux2", "Rorb", "Bcl11b", "Tbr1"],
    "Hippocampus": ["Prox1", "Calb1", "Wfs1", "Pcp4", "Grik4", "Bok"],
    "Striatum":    ["Drd1", "Drd2", "Penk", "Pdyn"],
}


def decode(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def build_gene_index(h5):
    """Return {gene_symbol: feature_slice_key_str} for genes with expression data."""
    names     = [decode(n) for n in h5["features"]["name"][:]]
    available = set(h5["feature_slices"].keys())
    return {name: str(i + 1) for i, name in enumerate(names)
            if str(i + 1) in available}, names


def read_gene_slice(h5, key):
    """Read /feature_slices/<key> and return DataFrame [x, y, umi]."""
    base = "/feature_slices/" + key
    rows = h5[base + "/row"][:]
    cols = h5[base + "/col"][:]
    data = h5[base + "/data"][:]
    return pd.DataFrame({"x": cols.astype(int),
                         "y": rows.astype(int),
                         "umi": data.astype(int)})


def plot_gene(df, gene, region, outdir, vmax_pct=99.5):
    total_umi = int(df["umi"].sum())
    n_nonzero = int((df["umi"] > 0).sum())
    print(f"  {gene}: total UMIs={total_umi:,}, non-zero spots={n_nonzero:,}")

    fig, ax = plt.subplots(figsize=(6, 6))

    if total_umi == 0:
        ax.text(0.5, 0.5, "no counts detected", transform=ax.transAxes,
                ha="center", va="center", color="grey", fontsize=11)
    else:
        vmax = float(np.percentile(df["umi"][df["umi"] > 0], vmax_pct))
        vmax = max(vmax, 1.0)

        n_pts = len(df)
        marker_area = min((6 * 6 * 72 ** 2) / max(n_pts, 1) * 0.5, 200)

        zero = df[df["umi"] == 0]
        pos  = df[df["umi"] > 0]

        if len(zero) > 0:
            ax.scatter(zero["x"], zero["y"], c="#d8d8d8", s=marker_area,
                       linewidths=0, rasterized=True, zorder=1)

        sc = ax.scatter(pos["x"], pos["y"], c=pos["umi"], s=marker_area,
                        cmap="Reds", vmin=0.5, vmax=vmax,
                        linewidths=0, rasterized=True, zorder=2)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("UMI count", fontsize=9)

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(gene + "  [" + region + "]  - real Visium HD",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("col (x)")
    ax.set_ylabel("row (y)")
    ax.tick_params(labelsize=8)
    fig.tight_layout()

    out_file = outdir / (region + "_" + gene + ".png")
    fig.savefig(str(out_file), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_file}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5",       required=True, type=Path,
                        help="Path to Visium HD feature_slice.h5")
    parser.add_argument("--outdir",   default=Path("real_visiumhd_heatmaps"), type=Path,
                        help="Output directory (default: real_visiumhd_heatmaps/)")
    parser.add_argument("--vmax-pct", default=99.5, type=float,
                        help="Percentile for colour-scale upper bound (default: 99.5)")
    args = parser.parse_args()

    if not args.h5.exists():
        sys.exit("ERROR: H5 file not found: " + str(args.h5))

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Opening " + str(args.h5) + " ...")
    with h5py.File(str(args.h5), "r") as h5:
        gene_index, all_names = build_gene_index(h5)
        print(f"  {len(all_names):,} features total")
        print(f"  {len(h5['feature_slices']):,} genes with non-zero expression")

        absent = [g for gs in MARKERS.values() for g in gs if g not in gene_index]
        if absent:
            print(f"  WARNING: genes not found or have no expression: {absent}")

        for region, genes in MARKERS.items():
            region_dir = args.outdir / region
            region_dir.mkdir(exist_ok=True)
            print(f"\n--- {region} ---")
            for gene in genes:
                if gene not in gene_index:
                    print(f"  SKIP {gene} (not in H5)")
                    continue
                df = read_gene_slice(h5, gene_index[gene])
                plot_gene(df, gene, region, region_dir, args.vmax_pct)

    print(f"\nDone. Figures saved to: {args.outdir}/")


if __name__ == "__main__":
    main()
