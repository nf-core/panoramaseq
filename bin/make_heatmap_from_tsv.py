#!/usr/bin/env python3
"""
Generate a spatial UMI heatmap from a starsolo_umi_heatmap.py data TSV.

Aggregates UMI counts at each unique (x,y) position, maps them onto a
regular 2D grid (preserving actual coordinate positions), and renders
with matplotlib imshow. This correctly handles simulation data where
multiple barcodes share the same spatial coordinate.

Usage:
  python3 make_heatmap_from_tsv.py \\
      --data  <heatmap_data.tsv> \\
      --out   <output.png> \\
      [--metric total_counts] [--percentile 99.5] [--title "My Sample"]
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm


def make_heatmap(data_tsv, out_png, title="UMI Heatmap", metric="total_counts",
                 colormap="viridis", percentile=99.5):

    print(f"Loading {data_tsv} ...")
    df = pd.read_csv(data_tsv, sep="\t")
    print(f"  {len(df):,} barcodes")

    # ── 1. Aggregate at unique (x,y) ──────────────────────────────────────────
    agg = df.groupby(["x", "y"])[metric].sum().reset_index()
    agg.columns = ["x", "y", "value"]
    n_unique = len(agg)
    print(f"  {n_unique:,} unique (x,y) positions  "
          f"({len(df)/n_unique:.1f} barcodes per position)")

    xs = agg["x"].values
    ys = agg["y"].values
    vals = agg["value"].values.astype(float)

    # ── 2. Colour normalisation ────────────────────────────────────────────────
    valid = vals[vals > 0]
    vmax = np.percentile(valid, percentile) if len(valid) > 0 else 1.0
    vmin = 0.0
    print(f"  Colour range: 0 – {vmax:.0f}  (p{percentile:.0f}={vmax:.0f})")

    # ── 3. Render with matplotlib scatter ─────────────────────────────────────
    # Compute spot size so each spot fills its territorial area.
    # Area per spot = canvas_area / n_spots → radius = sqrt(area/π)
    x_range = float(xs.max() - xs.min()) if len(xs) > 1 else 1.0
    y_range = float(ys.max() - ys.min()) if len(ys) > 1 else 1.0
    area_per_spot = (x_range * y_range) / n_unique

    dpi = 150
    fig_w, fig_h = 8, 8
    # Convert spot radius from data coords to matplotlib points²
    # One data unit in inches = fig_w / x_range, then × dpi = pixels per data unit
    px_per_data_unit_x = (fig_w * dpi) / (x_range * 1.05)
    px_per_data_unit_y = (fig_h * dpi) / (y_range * 1.05)
    px_per_data_unit = min(px_per_data_unit_x, px_per_data_unit_y)

    # Radius in data units that fills the territory; scale to matplotlib points²
    radius_data = np.sqrt(area_per_spot / np.pi)
    radius_pts  = radius_data * px_per_data_unit * (72.0 / dpi)  # pts
    marker_area = np.pi * radius_pts**2   # matplotlib s = area in points²
    print(f"  Marker radius: {radius_data:.1f} data units  →  {radius_pts:.1f} pts")

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = cm.get_cmap(colormap)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi, facecolor="black")
    ax.set_facecolor("black")

    sc = ax.scatter(
        xs, ys,
        c=vals,
        s=marker_area,
        marker="s",           # square markers
        cmap=cmap_obj,
        norm=norm,
        linewidths=0,
        rasterized=True,
    )

    cbar = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("UMI count", fontsize=8, color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    ax.set_xlim(xs.min() - radius_data * 2, xs.max() + radius_data * 2)
    ax.set_ylim(ys.min() - radius_data * 2, ys.max() + radius_data * 2)
    ax.set_title(title, fontsize=10, pad=6, color="white")
    ax.set_xlabel("x coordinate", fontsize=7, color="white")
    ax.set_ylabel("y coordinate", fontsize=7, color="white")
    ax.tick_params(colors="white", labelsize=6)

    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="black")
    plt.close()
    print(f"Saved: {out_png}")

    # ── 5. Stats ───────────────────────────────────────────────────────────────
    print(f"\n--- Summary ({metric}) ---")
    print(f"  Total barcodes:          {len(df):>10,}")
    print(f"  Unique spatial positions:{n_unique:>10,}")
    print(f"  Total UMIs (aggregated): {int(vals.sum()):>10,}")
    print(f"  Median UMIs per position:{int(np.median(vals[vals>0])) if len(valid)>0 else 0:>10,}")
    print(f"  Max UMIs per position:   {int(vals.max()) if len(vals)>0 else 0:>10,}")
    print(f"  Positions with UMI > 0:  {int((vals>0).sum()):>10,}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--data",       required=True,
                   help="heatmap data TSV from starsolo_umi_heatmap.py")
    p.add_argument("--out",        required=True, help="output PNG path")
    p.add_argument("--title",      default="UMI Heatmap")
    p.add_argument("--metric",     default="total_counts",
                   help="column to aggregate and plot (default: total_counts)")
    p.add_argument("--colormap",   default="viridis")
    p.add_argument("--percentile", type=float, default=99.5,
                   help="colour-scale clipping percentile (default 99.5)")
    args = p.parse_args()

    make_heatmap(
        data_tsv   = args.data,
        out_png    = args.out,
        title      = args.title,
        metric     = args.metric,
        colormap   = args.colormap,
        percentile = args.percentile,
    )


if __name__ == "__main__":
    main()
