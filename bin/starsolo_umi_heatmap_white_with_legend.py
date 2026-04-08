#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import scipy.io as sio
import anndata as ad
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import json
import gzip
import matplotlib.cm as mcm
import matplotlib.colors as mcolors

def load_starsolo_to_anndata(matrix_path, barcodes_path, features_path):
    """
    Load STARsolo output into AnnData object
    
    Args:
        matrix_path: Path to matrix.mtx.gz (gene x barcode count matrix)
        barcodes_path: Path to barcodes.tsv.gz (synthetic barcodes)
        features_path: Path to features.tsv.gz (gene annotations)
    
    Returns:
        AnnData object with synthetic barcodes as obs_names
    """
    print(f"Loading STARsolo data...")
    
    # Read matrix (Market Matrix format)
    matrix = sio.mmread(matrix_path).T.tocsr()
    
    # Read barcodes (synthetic 25bp)
    with gzip.open(barcodes_path, 'rt') as f:
        barcodes = pd.DataFrame({'barcode': [line.strip() for line in f]})
    
    # Read features
    with gzip.open(features_path, 'rt') as f:
        features = pd.read_csv(f, sep='\t', header=None)
        if features.shape[1] >= 2:
            features.columns = ['gene_id', 'gene_name'] + [f'col{i}' for i in range(2, features.shape[1])]
        else:
            features.columns = ['gene_id']
            features['gene_name'] = features['gene_id']
    
    # Create AnnData object
    adata = ad.AnnData(X=matrix, obs=barcodes, var=features[['gene_id', 'gene_name']])
    adata.var_names = adata.var['gene_id']
    adata.obs_names = adata.obs['barcode']
    
    # Calculate QC metrics
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes_by_counts'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata.var['total_counts'] = np.array(adata.X.sum(axis=0)).flatten()
    
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Total UMIs: {adata.X.sum():.0f}")
    
    return adata

def reverse_map_barcodes(adata, mapping_path):
    """
    Map synthetic barcodes back to original 36bp barcodes
    
    Args:
        adata: AnnData object with synthetic barcodes
        mapping_path: Path to barcode_mapping.tsv
    
    Returns:
        AnnData with original_barcode column added to obs
    """
    print(f"Reverse-mapping synthetic → original barcodes...")
    
    mapping = pd.read_csv(mapping_path, sep='\t')
    reverse_map = dict(zip(mapping['synthetic_barcode'], mapping['original_barcode']))
    
    adata.obs['original_barcode'] = adata.obs_names.map(reverse_map)
    
    # Filter out unmapped barcodes
    n_before = adata.n_obs
    adata = adata[~adata.obs['original_barcode'].isna()].copy()
    n_after = adata.n_obs
    
    print(f"  Mapped {n_after}/{n_before} barcodes to original sequences")
    
    return adata


def add_spatial_coordinates_direct(adata, coords_path):
    """
    Add spatial coordinates to AnnData object using SYNTHETIC barcodes (no reverse mapping needed).
    This is the optimized version that skips the expensive reverse mapping step.
    
    Args:
        adata: AnnData object with synthetic barcodes as obs_names
        coords_path: Path to barcode coordinates CSV with synthetic barcodes (cell, x, y)
    
    Returns:
        AnnData with spatial coordinates in obsm['spatial'], expanded to include all spots
    """
    print(f"Adding spatial coordinates (direct synthetic barcode matching)...")
    
    coords = pd.read_csv(coords_path)
    
    # Handle different column names
    if 'cell' in coords.columns:
        barcode_col = 'cell'
    elif 'barcode' in coords.columns:
        barcode_col = 'barcode'
    else:
        # Assume first column is barcode
        barcode_col = coords.columns[0]
    
    print(f"  Total spots in coordinate file: {len(coords)}")
    
    # Create a mapping from synthetic_barcode to counts
    barcode_to_data = {}
    for idx, row in adata.obs.iterrows():
        synthetic_bc = idx  # obs_names are already synthetic barcodes
        barcode_to_data[synthetic_bc] = {
            'total_counts': row['total_counts'],
            'n_genes_by_counts': row['n_genes_by_counts']
        }
    
    # Create new obs dataframe with ALL spots from coordinates
    new_obs = []
    spatial_coords = []
    
    for _, row in coords.iterrows():
        synthetic_bc = row[barcode_col]
        spatial_coords.append([row['x'], row['y']])
        
        if synthetic_bc in barcode_to_data:
            # Spot has UMI counts
            new_obs.append({
                'barcode': synthetic_bc,
                'total_counts': barcode_to_data[synthetic_bc]['total_counts'],
                'n_genes_by_counts': barcode_to_data[synthetic_bc]['n_genes_by_counts']
            })
        else:
            # Empty spot - add with zero counts
            new_obs.append({
                'barcode': synthetic_bc,
                'total_counts': 0,
                'n_genes_by_counts': 0
            })
    
    new_obs_df = pd.DataFrame(new_obs)
    new_obs_df.index = new_obs_df['barcode']
    
    # Create new AnnData with all spots
    # Initialize with zeros (we only care about obs metrics for heatmap)
    new_X = np.zeros((len(new_obs_df), adata.n_vars))
    
    # Fill in counts for spots that have data
    for i, barcode in enumerate(new_obs_df['barcode']):
        if barcode in adata.obs_names:
            # Find the row in original adata
            row_idx = adata.obs_names.get_loc(barcode)
            new_X[i, :] = adata.X[row_idx, :].toarray().flatten()
    
    # Create new AnnData
    new_adata = ad.AnnData(X=new_X, obs=new_obs_df, var=adata.var)
    new_adata.obsm['spatial'] = np.array(spatial_coords)
    
    spots_with_counts = (new_adata.obs['total_counts'] > 0).sum()
    print(f"  Total spots plotted: {new_adata.n_obs}")
    print(f"  Spots with UMI counts: {spots_with_counts}")
    print(f"  Empty spots: {new_adata.n_obs - spots_with_counts}")
    
    return new_adata


def add_spatial_coordinates(adata, coords_path):
    """
    Add spatial coordinates to AnnData object and expand to include ALL spots from coordinate file
    
    Args:
        adata: AnnData object with original_barcode in obs
        coords_path: Path to barcode coordinates CSV (cell, x, y)
    
    Returns:
        AnnData with spatial coordinates in obsm['spatial'], expanded to include all spots
    """
    print(f"Adding spatial coordinates...")
    
    coords = pd.read_csv(coords_path)
    
    # Handle different column names
    if 'cell' in coords.columns:
        coords.rename(columns={'cell': 'barcode'}, inplace=True)
    elif 'barcode' not in coords.columns:
        # Assume first column is barcode
        coords.columns = ['barcode'] + list(coords.columns[1:])
    
    print(f"  Total spots in coordinate file: {len(coords)}")
    
    # Create a mapping from original_barcode to counts
    barcode_to_data = {}
    for idx, row in adata.obs.iterrows():
        original_bc = row['original_barcode']
        barcode_to_data[original_bc] = {
            'total_counts': row['total_counts'],
            'n_genes_by_counts': row['n_genes_by_counts']
        }
    
    # Create new obs dataframe with ALL spots from coordinates
    new_obs = []
    spatial_coords = []
    
    for _, row in coords.iterrows():
        barcode = row['barcode']
        spatial_coords.append([row['x'], row['y']])
        
        if barcode in barcode_to_data:
            # Spot has UMI counts
            new_obs.append({
                'barcode': barcode,
                'original_barcode': barcode,
                'total_counts': barcode_to_data[barcode]['total_counts'],
                'n_genes_by_counts': barcode_to_data[barcode]['n_genes_by_counts']
            })
        else:
            # Empty spot - add with zero counts
            new_obs.append({
                'barcode': barcode,
                'original_barcode': barcode,
                'total_counts': 0,
                'n_genes_by_counts': 0
            })
    
    new_obs_df = pd.DataFrame(new_obs)
    new_obs_df.index = new_obs_df['barcode']
    
    # Create new AnnData with all spots
    # Initialize with zeros (we only care about obs metrics for heatmap)
    new_X = np.zeros((len(new_obs_df), adata.n_vars))
    
    # Fill in counts for spots that have data
    for i, barcode in enumerate(new_obs_df['barcode']):
        if barcode in adata.obs['original_barcode'].values:
            # Find the row in original adata
            orig_idx = adata.obs[adata.obs['original_barcode'] == barcode].index[0]
            row_idx = adata.obs.index.get_loc(orig_idx)
            new_X[i, :] = adata.X[row_idx, :].toarray().flatten()
    
    # Create new AnnData
    new_adata = ad.AnnData(X=new_X, obs=new_obs_df, var=adata.var)
    new_adata.obsm['spatial'] = np.array(spatial_coords)
    
    spots_with_counts = (new_adata.obs['total_counts'] > 0).sum()
    print(f"  Total spots plotted: {new_adata.n_obs}")
    print(f"  Spots with UMI counts: {spots_with_counts}")
    print(f"  Empty spots: {new_adata.n_obs - spots_with_counts}")
    
    return new_adata

def bin_spatial_data(adata, bin_size=1):
    """
    Aggregate spots into spatial bins
    
    Args:
        adata: AnnData with spatial coordinates in obsm['spatial']
        bin_size: Number of spots per bin dimension (1=no binning, 10=10x10 bins)
    
    Returns:
        New AnnData with binned data
    """
    if bin_size <= 1:
        # No binning, just add metadata
        adata.uns['bin_size'] = 1
        adata.uns['spots_per_bin'] = 1
        return adata
    
    print(f"Binning spots into {bin_size}x{bin_size} bins...")
    
    spatial = adata.obsm['spatial']
    x_coords = spatial[:, 0]
    y_coords = spatial[:, 1]
    
    # Calculate bin boundaries
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    # Assign each spot to a bin using floor division
    # This ensures each spot belongs to exactly ONE bin
    x_bin_idx = ((x_coords - x_min) // bin_size).astype(int)
    y_bin_idx = ((y_coords - y_min) // bin_size).astype(int)
    
    adata.obs['x_bin'] = x_bin_idx
    adata.obs['y_bin'] = y_bin_idx
    adata.obs['bin_id'] = adata.obs['x_bin'].astype(str) + '_' + adata.obs['y_bin'].astype(str)
    
    # Group by bin and aggregate
    bin_ids = adata.obs['bin_id'].unique()
    
    binned_matrices = []
    binned_obs = []
    binned_spatial = []
    
    for bin_id in bin_ids:
        mask = adata.obs['bin_id'] == bin_id
        bin_data = adata[mask]
        
        # Sum counts across spots in bin
        summed_counts = bin_data.X.sum(axis=0)
        binned_matrices.append(summed_counts)
        
        # Use bin indices to calculate GRID-ALIGNED center coordinates
        # This prevents overlapping by ensuring bins are on a regular grid
        x_bin = bin_data.obs['x_bin'].iloc[0]
        y_bin = bin_data.obs['y_bin'].iloc[0]
        
        # Calculate bin center: min_coord + (bin_index * bin_size) + (bin_size / 2)
        bin_x = x_min + (x_bin * bin_size) + (bin_size / 2.0)
        bin_y = y_min + (y_bin * bin_size) + (bin_size / 2.0)
        binned_spatial.append([bin_x, bin_y])
        
        # Store metadata
        binned_obs.append({
            'bin_id': bin_id,
            'total_counts': summed_counts.sum(),
            'n_genes_by_counts': (summed_counts > 0).sum(),
            'n_spots': mask.sum(),
            'x_bin': x_bin,
            'y_bin': y_bin,
            'bin_x_center': bin_x,
            'bin_y_center': bin_y
        })
    
    # Create new AnnData with binned data
    binned_X = np.vstack(binned_matrices)
    binned_obs_df = pd.DataFrame(binned_obs)
    binned_obs_df.index = binned_obs_df['bin_id']
    
    adata_binned = ad.AnnData(X=binned_X, obs=binned_obs_df, var=adata.var)
    adata_binned.obsm['spatial'] = np.array(binned_spatial)
    adata_binned.uns['bin_size'] = bin_size
    adata_binned.uns['spots_per_bin'] = bin_size
    
    # Validate no overlapping bins (each bin should have unique coordinates)
    spatial_coords = adata_binned.obsm['spatial']
    unique_coords = np.unique(spatial_coords, axis=0)
    if len(unique_coords) < len(spatial_coords):
        print(f"  WARNING: Found {len(spatial_coords) - len(unique_coords)} duplicate bin coordinates!")
        print(f"  This may indicate overlapping bins in visualization.")
    else:
        print(f"  Validation: All {len(spatial_coords)} bins have unique coordinates ✓")
    
    print(f"  Created {adata_binned.n_obs} bins from {adata.n_obs} spots")
    print(f"  Bins are grid-aligned to prevent overlap")
    
    return adata_binned

def square_heatmap(
    adata, 
    metric='total_counts',
    percentile=99.5,
    filename='square_heatmap.png',
    legend_filename=None,
    spot_size=None,
    colormap='viridis',
    title=None
):
    """
    Generate square-spot heatmap using PIL (adapted from hexbin_plot)
    
    Args:
        adata: AnnData with spatial coordinates and metric in obs
        metric: Column name in adata.obs to visualize (or 'rgb' for 3-channel)
        percentile: Percentile for color clipping (default: 99.5)
        filename: Output filename
        legend_filename: Optional output path for a legend PNG
        spot_size: Size of each square spot in pixels (auto if None)
        colormap: Color scheme ('viridis', 'magma', 'plasma', 'rgb')
        title: Plot title
    """
    print(f"Generating square heatmap...")
    
    cap = None
    if metric == 'rgb':
        # RGB mode: use total_counts (R), n_genes_by_counts (G), blue=0
        r_scores = adata.obs['total_counts'].values
        g_scores = adata.obs['n_genes_by_counts'].values
        b_scores = np.zeros(len(r_scores))
        title = title or 'RGB Heatmap (R=UMI, G=Genes, B=0)'
    else:
        # Single metric mode
        scores = adata.obs[metric].values
        title = title or f'Heatmap of {metric}'
    
    spatial = adata.obsm['spatial'].astype(int)
    
    # Determine spot size
    if spot_size is None:
        # Auto-calculate based on bin size
        bin_size = adata.uns.get('bin_size', 1)
        # For binned data, use bin_size as the spot_size to fill the bin area
        # For unbinned data, use a small fixed size
        if bin_size > 1:
            spot_size = int(bin_size)  # Each bin fills its allocated area
        else:
            spot_size = 2  # Small spots for unbinned data
    
    print(f"  Spot/bin size: {spot_size} pixels (bin_size={adata.uns.get('bin_size', 1)})")
    
    # Calculate image dimensions
    max_x, max_y = spatial.max(axis=0)
    min_x, min_y = spatial.min(axis=0)
    
    img_width = max_x - min_x + spot_size * 2
    img_height = max_y - min_y + spot_size * 2
    
    print(f"  Image size: {img_width}x{img_height} pixels")
    print(f"  Spot size: {spot_size} pixels")
    print(f"  Number of spots: {len(spatial)}")
    
    final_img = Image.new('RGB', (img_width, img_height), (255, 255, 255))
    draw = ImageDraw.Draw(final_img)
    
    if metric == 'rgb':
        # RGB mode
        r_norm = np.clip((r_scores / np.percentile(r_scores, percentile) * 255).astype(int), 0, 255)
        g_norm = np.clip((g_scores / np.percentile(g_scores, percentile) * 255).astype(int), 0, 255)
        b_norm = np.clip((b_scores / np.percentile(b_scores, percentile) * 255).astype(int), 0, 255)
        
        for coord, r, g, b in zip(spatial, r_norm, g_norm, b_norm):
            x = int(coord[0] - min_x)
            y = int(max_y - coord[1])  # Invert Y axis
            
            # Draw square
            draw.rectangle(
                [x - spot_size//2, y - spot_size//2, x + spot_size//2, y + spot_size//2],
                fill=(int(r), int(g), int(b))
            )
    else:
        # Single channel with colormap
        # Handle edge case: all zeros or very low values
        # Log1p scaling: compresses the high-count tail and lifts low counts
        # away from zero, so every expressed spot gets a clearly visible colour.
        scores_log = np.log1p(scores.astype(float))
        cap = np.log1p(np.percentile(scores[scores > 0], percentile)
                   if (scores > 0).any() else 1.0)
        if cap > 0:
            scores_norm = scores_log / cap
        else:
            scores_norm = np.zeros_like(scores_log)

        scores_norm = np.clip(scores_norm, 0, 1)

        # Apply colormap
        colors = apply_colormap(scores_norm, colormap)
        
        for coord, color in zip(spatial, colors):
            if color is None:
                continue      # zero-count spot — leave white background
            x = int(coord[0] - min_x)
            y = int(max_y - coord[1])  # Invert Y axis

            draw.rectangle(
                [x - spot_size//2, y - spot_size//2, x + spot_size//2, y + spot_size//2],
                fill=color
            )
    
    final_img.save(filename)
    print(f"  Saved heatmap: {filename}")

    if legend_filename:
        if metric == 'rgb':
            print("  Skipping legend: RGB mode does not use a single scalar color scale")
        else:
            save_color_legend(
                out_path=legend_filename,
                colormap=colormap,
                percentile=percentile,
                cap_log=float(cap) if cap is not None else 0.0,
                metric=metric,
            )
            print(f"  Saved legend: {legend_filename}")
    
    return final_img

def apply_colormap(values, colormap='Reds'):
    """
    Apply a matplotlib colormap to normalised values [0, 1].
    Zero values are returned as None so the caller can skip drawing them
    (leaving the white background visible).

    Returns:
        List of (R, G, B) uint8 tuples, or None for zero-count spots.
    """
    values = np.nan_to_num(values, nan=0.0)

    # Remap non-zero values to [0.25, 1.0] so no expressed spot maps to the
    # near-white end of the colormap. Zero stays at 0 (→ white background).
    display = np.where(values > 0, 0.25 + 0.75 * values, 0.0)

    cmap = mcm.get_cmap(colormap)
    rgba = cmap(display)                       # (N, 4) float 0-1
    rgb  = (rgba[:, :3] * 255).astype(np.uint8)

    result = []
    for v, c in zip(values, rgb):
        result.append(None if v == 0.0 else tuple(c))
    return result


def save_color_legend(out_path, colormap, percentile, cap_log, metric):
    """
    Save a standalone vertical legend that matches the heatmap colour mapping.
    """
    bar_h = 600
    bar_w = 64
    pad = 24
    text_x = pad + bar_w + 18
    img_w = 360
    img_h = bar_h + pad * 2

    img = Image.new('RGB', (img_w, img_h), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    y_norm = np.linspace(1.0, 0.0, bar_h)
    bar_colors = apply_colormap(y_norm, colormap)
    for i, col in enumerate(bar_colors):
        color = (255, 255, 255) if col is None else col
        y0 = pad + i
        draw.rectangle([pad, y0, pad + bar_w, y0], fill=color)

    draw.rectangle([pad, pad, pad + bar_w, pad + bar_h], outline=(0, 0, 0), width=1)

    draw.text((text_x, pad - 14), f"{metric} (log1p scale)", fill=(0, 0, 0), font=font)
    draw.text((text_x, pad + 4), f"clip at p{percentile:g}", fill=(0, 0, 0), font=font)

    tick_fracs = [0.0, 0.25, 0.5, 0.75, 1.0]
    for f in tick_fracs:
        y = int(pad + bar_h * (1.0 - f))
        draw.line([pad + bar_w + 3, y, pad + bar_w + 12, y], fill=(0, 0, 0), width=1)
        if f == 0.0:
            label = '0'
        else:
            label = f"{int(round(np.expm1(f * max(cap_log, 0.0)))):,}"
        draw.text((text_x, y - 7), label, fill=(0, 0, 0), font=font)

    suffix = Path(out_path).suffix.lower()
    if suffix == '.pdf':
        img.save(out_path, 'PDF', resolution=300.0)
    else:
        img.save(out_path)

def main():
    parser = argparse.ArgumentParser(
        description='Generate spatial UMI count heatmap from STARsolo output'
    )
    parser.add_argument('--matrix', required=True, help='Path to matrix.mtx.gz')
    parser.add_argument('--barcodes', required=True, help='Path to barcodes.tsv.gz')
    parser.add_argument('--features', required=True, help='Path to features.tsv.gz')
    parser.add_argument('--mapping', required=False, help='Path to barcode_mapping.tsv (optional if coords has synthetic barcodes)')
    parser.add_argument('--coords', required=True, help='Path to spatial coordinates CSV')
    parser.add_argument('--output-png', required=True, help='Output PNG heatmap')
    parser.add_argument('--output-legend',
                       help='Optional output path for color legend (e.g. .png or .pdf)')
    parser.add_argument('--output-h5ad', help='Optional: save AnnData object')
    parser.add_argument('--output-data', required=True, help='Output data TSV')
    parser.add_argument('--output-stats', required=True, help='Output stats JSON')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--bin-size', type=int, default=1,
                       help='Spatial bin size (1=no binning, 10=10x10 bins)')
    parser.add_argument('--metric', default='total_counts',
                       help='Metric to visualize (total_counts, n_genes_by_counts, or rgb)')
    parser.add_argument('--colormap', default='Reds',
                       choices=['Reds', 'Blues', 'Greens', 'Oranges', 'Purples',
                                'YlOrRd', 'YlGnBu', 'RdPu', 'hot_r', 'plasma',
                                'viridis', 'rgb'],
                       help='Matplotlib colormap (default: Reds — works well on white)')
    parser.add_argument('--percentile', type=float, default=99.5,
                       help='Percentile for color clipping')
    parser.add_argument('--spot-size', type=int, default=None,
                       help='Size of each spot in pixels (auto if not specified)')
    
    args = parser.parse_args()
    
    print(f"=== STARsolo UMI Heatmap Generator ===")
    print(f"Sample: {args.sample_id}")
    print(f"Bin size: {args.bin_size}")
    print(f"Metric: {args.metric}")
    print(f"Colormap: {args.colormap}")
    print(f"Mode: {'Direct synthetic coords' if not args.mapping else 'With reverse mapping'}")
    print(f"==========================================")
    
    # 1. Load STARsolo data into AnnData
    adata = load_starsolo_to_anndata(args.matrix, args.barcodes, args.features)
    
    # 2. Add spatial coordinates
    if args.mapping:
        # OLD MODE: Reverse-map synthetic → original barcodes, then match with coords
        print("Using legacy mode: reverse mapping synthetic barcodes...")
        adata = reverse_map_barcodes(adata, args.mapping)
        adata = add_spatial_coordinates(adata, args.coords)
    else:
        # NEW MODE: Direct matching with synthetic coords (MUCH FASTER)
        print("Using optimized mode: direct synthetic barcode matching...")
        adata = add_spatial_coordinates_direct(adata, args.coords)
    
    # 4. Spatial binning (if requested)
    if args.bin_size > 1:
        adata = bin_spatial_data(adata, args.bin_size)
    else:
        adata.uns['bin_size'] = 1
    
    # 5. Generate heatmap
    metric = args.metric if args.colormap != 'rgb' else 'rgb'
    square_heatmap(
        adata,
        metric=metric,
        percentile=args.percentile,
        filename=args.output_png,
        legend_filename=args.output_legend,
        spot_size=args.spot_size,
        colormap=args.colormap,
        title=f"{args.sample_id} - UMI Count Heatmap"
    )
    
    # 6. Save data table
    output_df = pd.DataFrame({
        'spot_id': adata.obs_names,
        'x': adata.obsm['spatial'][:, 0],
        'y': adata.obsm['spatial'][:, 1],
        'total_counts': adata.obs['total_counts'],
        'n_genes': adata.obs['n_genes_by_counts']
    })
    
    if args.bin_size > 1:
        output_df['n_spots_in_bin'] = adata.obs['n_spots']
    
    output_df.to_csv(args.output_data, sep='\t', index=False)
    print(f"Saved data table: {args.output_data}")
    
    # 7. Calculate statistics
    stats = {
        'sample_id': args.sample_id,
        'bin_size': args.bin_size,
        'n_spots': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'total_umis': int(adata.X.sum()),
        'mean_umi_per_spot': float(adata.obs['total_counts'].mean()),
        'median_umi_per_spot': float(adata.obs['total_counts'].median()),
        'max_umi': int(adata.obs['total_counts'].max()),
        'min_umi': int(adata.obs['total_counts'].min()),
        'mean_genes_per_spot': float(adata.obs['n_genes_by_counts'].mean())
    }
    
    with open(args.output_stats, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Saved statistics: {args.output_stats}")
    
    # 8. Optional: save AnnData object
    if args.output_h5ad:
        adata.write_h5ad(args.output_h5ad, compression='gzip')
        print(f"Saved AnnData: {args.output_h5ad}")
    
    print("=== Heatmap Generation Complete ===")
    print(f"Total spots: {stats['n_spots']}")
    print(f"Total UMIs: {stats['total_umis']:,}")
    print(f"Mean UMIs/spot: {stats['mean_umi_per_spot']:.1f}")
    print("===================================")

if __name__ == '__main__':
    main()
