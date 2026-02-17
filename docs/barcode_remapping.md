# Barcode Remapping for STARsolo

## Overview

STARsolo has a hard-coded limit of 31bp for cell barcodes. Since panoramaseq uses 36bp barcodes, we implement a **synthetic barcode remapping** strategy to bypass this limitation without losing barcode information.

## Approach

1. **QUIK** generates a whitelist of corrected 36bp barcodes
2. **REMAP_BARCODES** creates unique synthetic ≤31bp barcodes (default: 25bp) for each original barcode
3. **FASTQ files** are rewritten with synthetic barcodes replacing original ones
4. **STARsolo** processes data using synthetic barcodes
5. **Mapping table** is preserved for downstream restoration

## Why 25bp synthetic barcodes?

- Well below STARsolo's 31bp limit (safe margin)
- Can encode up to 4^25 unique barcodes (>10^15), far exceeding typical spatial assay capacity
- Smaller barcodes = faster alignment
- Deterministic encoding (same whitelist → same synthetic barcodes)

## Implementation

### Synthetic Barcode Generation

Barcodes are generated using base-4 encoding:
- Each barcode gets a sequential index (0, 1, 2, ...)
- Index is converted to base-4 (DNA alphabet: A=0, C=1, G=2, T=3)
- Result is left-padded with 'A' to reach target length (25bp)

Example:
```
Index 0   → AAAAAAAAAAAAAAAAAAAAAAAAA (25 A's)
Index 1   → AAAAAAAAAAAAAAAAAAAAAAAC
Index 2   → AAAAAAAAAAAAAAAAAAAAAAAG
Index 3   → AAAAAAAAAAAAAAAAAAAAAAAT
Index 4   → AAAAAAAAAAAAAAAAAAAAAACA
...
```

### File Outputs

1. **`*_synthetic.fastq.gz`** - FASTQ with synthetic barcodes
2. **`*_whitelist_synthetic.txt`** - Whitelist of synthetic barcodes
3. **`*_barcode_mapping.tsv`** - Mapping table (original ↔ synthetic)

### Mapping Table Format

```tsv
original_barcode	synthetic_barcode
ACGTACGTACGTACGTACGTACGTACGTACGTACGT	AAAAAAAAAAAAAAAAAAAAAAAAA
TGCATGCATGCATGCATGCATGCATGCATGCATGCA	AAAAAAAAAAAAAAAAAAAAAAAC
...
```

## Restoring Original Barcodes

### Option 1: Post-processing STARsolo Output

After STARsolo completes, use the mapping table to restore original barcodes in the output matrices:

```bash
# Restore barcodes in STARsolo output
python3 << 'EOF'
import pandas as pd
import gzip

# Load mapping
mapping = pd.read_csv('barcode_mapping.tsv', sep='\t', index_col='synthetic_barcode')
reverse_map = mapping['original_barcode'].to_dict()

# Read STARsolo barcodes.tsv
bc_file = 'Solo.out/Gene/filtered/barcodes.tsv.gz'
with gzip.open(bc_file, 'rt') as f:
    synthetic_bcs = [line.strip() for line in f]

# Restore original barcodes
original_bcs = [reverse_map.get(bc, bc) for bc in synthetic_bcs]

# Write restored barcodes
with gzip.open('barcodes_original.tsv.gz', 'wt') as f:
    for bc in original_bcs:
        f.write(f"{bc}\n")
EOF
```

### Option 2: Integration in STARSOLO_TO_H5AD

The `STARSOLO_TO_H5AD` module could be modified to automatically apply the reverse mapping:

```python
# In STARSOLO_TO_H5AD module
import pandas as pd

# Load mapping if provided
if mapping_file:
    mapping_df = pd.read_csv(mapping_file, sep='\t', index_col='synthetic_barcode')
    reverse_map = mapping_df['original_barcode'].to_dict()
    
    # Restore barcodes in AnnData object
    adata.obs.index = [reverse_map.get(bc, bc) for bc in adata.obs.index]
```

## Validation

The remapping process includes several safety checks:

- ✅ All whitelist barcodes have unique synthetic mappings
- ✅ No synthetic barcode collisions
- ✅ All FASTQ barcodes exist in whitelist
- ✅ Read counts match between original and remapped FASTQ
- ✅ Synthetic barcode length ≤ 31bp

## Performance

- **Memory**: Streaming FASTQ processing, minimal memory footprint
- **Speed**: ~100,000 reads/second (single-threaded Python)
- **Determinism**: Same whitelist → same synthetic barcodes (reproducible)

## Configuration

Synthetic barcode length can be adjusted in `conf/modules.config`:

```groovy
withName: '.*:REMAP_BARCODES_FOR_STARSOLO' {
    ext.synthetic_length = 25   // Default: 25bp (max: 31bp)
}
```

**Note**: Changing `synthetic_length` affects STARsolo parameters:
- `--soloCBlen` must match synthetic barcode length
- `--soloUMIstart` must be `synthetic_length + 1`

## Comparison with Truncation

We tested truncating 36bp barcodes to 31bp:

| Metric | Full 36bp | Truncated 31bp | Synthetic 25bp |
|--------|-----------|----------------|----------------|
| Unique barcodes | 14,684 | 13,343 (-9.1%) | 14,684 (✅) |
| Information loss | None | Last 5bp lost | None (reversible) |
| STARsolo compatible | ❌ No | ✅ Yes | ✅ Yes |
| Downstream restoration | N/A | ❌ Impossible | ✅ Full restoration |

**Conclusion**: Synthetic barcode remapping preserves all 36bp barcode information while maintaining STARsolo compatibility.
