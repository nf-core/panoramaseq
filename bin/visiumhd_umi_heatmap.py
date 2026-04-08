#!/usr/bin/env python3
"""
Generate a spatial UMI-count heatmap directly from a Visium HD feature_slice.h5 file,
using the same square-pixel rendering approach as starsolo_umi_heatmap.py.

The feature_slice.h5 stores every expressed gene as a sparse slice:
  /feature_slices/<1-based-idx>/row   – tissue row indices
  /feature_slices/<1-based-idx>/col   – tissue col indices
  /feature_slices/<1-based-idx>/data  – UMI counts

This script:
  1. Iterates all feature slices to accumulate total_umis and n_genes per spot.
  2. Renders a square-pixel image (one pixel = one 2-um Visium HD bin) using PIL,
     the same colormap logic as starsolo_umi_heatmap.py.
  3. Saves PNG, TSV data table, and JSON stats.

Usage:
    python3 visiumhd_umi_heatmap.py \\
        --h5      <path>/feature_slice.h5 \\
        --outdir  real_visiumhd_heatmap/ \\
        --sample  VisiumHD_MouseBrain \\
        [--bin-size 4] [--metric total_counts] [--colormap viridis]
"""

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

# ---------------------------------------------------------------------------
# Colormap (identical to starsolo_umi_heatmap.py so outputs are comparable)
# ---------------------------------------------------------------------------

def apply_colormap(values, colormap="viridis"):
    values = np.nan_to_num(values, nan=0.0)
    colors = np.zeros((len(values), 3), dtype=np.uint8)

    if colormap == "viridis":
        for i, v in enumerate(values):
            colors[i] = [
                int(255 * (0.267 + 0.3 * v + 0.43 * v ** 2)),
                int(255 * (0.004 + 0.5 * v + 0.496 * v ** 2)),
                int(255 * (0.329 + 0.7 * v - 1.029 * v ** 2)),
            ]
    elif colormap == "magma":
        for i, v in enumerate(values):
            colors[i] = [
                int(255 * v ** 0.5),
                int(255 * v ** 2),
                int(255 * v ** 1.5),
            ]
    elif colormap == "plasma":
        for i, v in enumerate(values):
            colors[i] = [
                int(255 * (0.5 + 0.5 * v)),
                int(255 * v ** 2),
                int(255 * (0.8 - 0.3 * v)),
            ]
    elif colormap == "hot":
        for i, v in enumerate(values):
            if v < 0.33:
                colors[i] = [int(255 * v * 3), 0, 0]
            elif v < 0.67:
                colors[i] = [255, int(255 * (v - 0.33) * 3), 0]
            else:
                colors[i] = [255, 255, int(255 * (v - 0.67) * 3)]
    else:
        gray = (values * 255).astype(np.uint8)
        colors = np.column_stack([gray, gray, gray])

    return np.clip(colors, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------------------
# H5 reading
# ---------------------------------------------------------------------------

def accumulate_spots(h5_path: str) -> pd.DataFrame:
    """
    Iterate all /feature_slices in the H5 and return a DataFrame:
        col (x), row (y), total_umis, n_genes

    Uses numpy bincount for O(N) vectorised accumulation — no Python inner loop.
    Strategy:
      1. First pass to find max row/col bounds (reads just one array per slice).
      2. Allocate two flat numpy arrays indexed by (row * stride + col).
      3. Second pass: for each slice, compute linear indices and call np.add.at.
    """
    with h5py.File(h5_path, "r") as h5:
        slices   = list(h5["feature_slices"].keys())
        n_slices = len(slices)
        print(f"  Scanning {n_slices:,} feature slices (pass 1: bounds) ...")

        max_row = max_col = 0
        for i, key in enumerate(slices):
            if (i + 1) % 4000 == 0:
                print(f"    bounds pass {i+1:,}/{n_slices:,} ...")
            base = f"/feature_slices/{key}"
            r = h5[f"{base}/row"][:]
            c = h5[f"{base}/col"][:]
            if len(r):
                if r.max() > max_row: max_row = int(r.max())
                if c.max() > max_col: max_col = int(c.max())

        stride = max_col + 1
        n_flat = (max_row + 1) * stride
        print(f"  Grid: rows 0-{max_row}, cols 0-{max_col} → {n_flat:,} flat cells")

        umi_flat  = np.zeros(n_flat, dtype=np.int64)
        gene_flat = np.zeros(n_flat, dtype=np.int32)

        print(f"  Pass 2: accumulating UMIs ...")
        for i, key in enumerate(slices):
            if (i + 1) % 4000 == 0:
                print(f"    accum pass {i+1:,}/{n_slices:,} ...")
            base = f"/feature_slices/{key}"
            r    = h5[f"{base}/row"][:].astype(np.int64)
            c    = h5[f"{base}/col"][:].astype(np.int64)
            d    = h5[f"{base}/data"][:].astype(np.int64)
            idx  = r * stride + c
            np.add.at(umi_flat,  idx, d)
            np.add.at(gene_flat, idx, 1)

    nonzero = np.flatnonzero(umi_flat)
    print(f"  Done. {len(nonzero):,} unique spots with ≥1 UMI.")

    rows_out = nonzero // stride
    cols_out = nonzero  % stride
    df = pd.DataFrame({
        "x":          cols_out.astype(int),
        "y":          rows_out.astype(int),
        "total_umis": umi_flat[nonzero],
        "n_genes":    gene_flat[nonzero],
    })
    return df


# ---------------------------------------------------------------------------
# Spatial binning (same logic as starsolo_umi_heatmap.py)
# ---------------------------------------------------------------------------

def bin_spots(df: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    if bin_size <= 1:
        return df.copy()

    print(f"  Binning into {bin_size}x{bin_size} pixel bins ...")
    x_min, y_min = df["x"].min(), df["y"].min()
    df = df.copy()
    df["xb"] = ((df["x"] - x_min) // bin_size).astype(int)
    df["yb"] = ((df["y"] - y_min) // bin_size).astype(int)

    agg = df.groupby(["xb", "yb"], as_index=False).agg(
        total_umis=("total_umis", "sum"),
        n_genes=("n_genes", "sum"),
        n_spots=("total_umis", "count"),
    )
    agg["x"] = x_min + agg["xb"] * bin_size + bin_size / 2
    agg["y"] = y_min + agg["yb"] * bin_size + bin_size / 2
    print(f"  {len(agg):,} bins from {len(df):,} spots.")
    return agg[["x", "y", "total_umis", "n_genes", "n_spots"]]


# ---------------------------------------------------------------------------
# Square-pixel heatmap (same approach as starsolo_umi_heatmap.py)
# ---------------------------------------------------------------------------

def square_heatmap(df: pd.DataFrame, metric: str, bin_size: int,
                   percentile: float, colormap: str,
                   out_png: Path, title: str = "") -> None:

    scores = df[metric].values.astype(float)
    x = df["x"].values.astype(int)
    y = df["y"].values.astype(int)

    spot_size = max(bin_size, 2)
    x_min, y_min = x.min(), y.min()
    x_max, y_max = x.max(), y.max()

    img_w = x_max - x_min + spot_size * 2
    img_h = y_max - y_min + spot_size * 2
    print(f"  Image: {img_w}x{img_h} px  spot_size={spot_size}px  n_spots={len(df):,}")

    pct_val = np.percentile(scores[scores > 0], percentile) if (scores > 0).any() else 1.0
    scores_norm = np.clip(scores / max(pct_val, 1.0), 0.0, 1.0)
    colors = apply_colormap(scores_norm, colormap)

    img = Image.new("RGB", (img_w, img_h), (0, 0, 0))
    draw = ImageDraw.Draw(img)
    half = spot_size // 2

    for xi, yi, col in zip(x, y, colors):
        px = int(xi - x_min)
        py = int(y_max - yi)         # invert Y (image convention)
        draw.rectangle(
            [px - half, py - half, px + half, py + half],
            fill=tuple(col),
        )

    img.save(str(out_png))
    print(f"  Saved: {out_png}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--h5",       required=True,  help="Path to feature_slice.h5")
    parser.add_argument("--outdir",   default="real_visiumhd_heatmap",
                        help="Output directory (default: real_visiumhd_heatmap/)")
    parser.add_argument("--sample",   default="VisiumHD",
                        help="Sample name used in filenames and title")
    parser.add_argument("--bin-size", type=int,   default=1,
                        help="Aggregate NxN Visium HD bins (1=no aggregation, 4=8um bins)")
    parser.add_argument("--metric",   default="total_umis",
                        choices=["total_umis", "n_genes"],
                        help="Metric to colour (default: total_umis)")
    parser.add_argument("--colormap", default="viridis",
                        choices=["viridis", "magma", "plasma", "hot"],
                        help="Colour scheme (default: viridis)")
    parser.add_argument("--percentile", type=float, default=99.5,
                        help="Colour-scale upper-bound percentile (default: 99.5)")
    args = parser.parse_args()

    h5_path = Path(args.h5)
    if not h5_path.exists():
        sys.exit(f"ERROR: {h5_path} not found")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"=== Visium HD UMI Heatmap ===")
    print(f"H5:       {h5_path}")
    print(f"Sample:   {args.sample}")
    print(f"Bin size: {args.bin_size}")
    print(f"Metric:   {args.metric}")
    print(f"Colormap: {args.colormap}")
    print(f"Outdir:   {outdir}")
    print()

    # 1. Accumulate per-spot counts from all feature slices
    print("Step 1: accumulating UMIs per spot ...")
    df = accumulate_spots(str(h5_path))

    # 2. Bin if requested
    print("Step 2: spatial binning ...")
    df = bin_spots(df, args.bin_size)

    # 3. Render heatmap
    print("Step 3: rendering heatmap ...")
    tag      = f"{args.sample}_bin{args.bin_size}_{args.metric}_{args.colormap}"
    out_png  = outdir / f"{tag}.png"
    square_heatmap(df, metric=args.metric, bin_size=args.bin_size,
                   percentile=args.percentile, colormap=args.colormap,
                   out_png=out_png,
                   title=f"{args.sample} – {args.metric} (bin={args.bin_size})")

    # 4. Save data table
    out_tsv = outdir / f"{tag}_data.tsv"
    df.to_csv(str(out_tsv), sep="\t", index=False)
    print(f"  Saved data table: {out_tsv}")

    # 5. Save stats JSON (mirrors starsolo_umi_heatmap.py output format)
    nonzero = df[df["total_umis"] > 0]
    stats = {
        "sample_id":           args.sample,
        "bin_size":            args.bin_size,
        "n_spots_total":       int(len(df)),
        "n_spots_with_umis":   int(len(nonzero)),
        "total_umis":          int(df["total_umis"].sum()),
        "mean_umi_per_spot":   float(df["total_umis"].mean()),
        "median_umi_per_spot": float(df["total_umis"].median()),
        "max_umi":             int(df["total_umis"].max()),
        "mean_genes_per_spot": float(df["n_genes"].mean()),
        "percentile_used":     args.percentile,
        "metric":              args.metric,
        "colormap":            args.colormap,
    }
    out_json = outdir / f"{tag}_stats.json"
    with open(str(out_json), "w") as fh:
        json.dump(stats, fh, indent=2)
    print(f"  Saved stats: {out_json}")

    print()
    print("=== Summary ===")
    print(f"Total spots (with UMIs): {stats['n_spots_with_umis']:,}")
    print(f"Total UMIs:              {stats['total_umis']:,}")
    print(f"Mean UMIs/spot:          {stats['mean_umi_per_spot']:.2f}")
    print(f"Median UMIs/spot:        {stats['median_umi_per_spot']:.1f}")
    print(f"Max UMIs/spot:           {stats['max_umi']:,}")
    print(f"Mean genes/spot:         {stats['mean_genes_per_spot']:.2f}")


if __name__ == "__main__":
    main()
