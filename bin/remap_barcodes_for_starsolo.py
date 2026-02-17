#!/usr/bin/env python3
"""
Remap 36bp barcodes to synthetic ≤31bp barcodes for STARsolo compatibility.

This script:
1. Reads a whitelist of 36bp barcodes (from QUIK)
2. Generates unique synthetic barcodes (≤31bp) for each
3. Creates a mapping table (original -> synthetic)
4. Rewrites FASTQ reads with synthetic barcodes
5. Rewrites whitelist with synthetic barcodes

Author: nf-core/panoramaseq
"""

import gzip
import sys
import argparse
from pathlib import Path


def encode_index_to_dna(index, length=25):
    """
    Encode an integer index as a DNA sequence of fixed length.
    Uses base-4 encoding: A=0, C=1, G=2, T=3.
    
    Args:
        index: Integer to encode
        length: Target DNA sequence length (default 25bp, max 31bp)
    
    Returns:
        DNA sequence string of specified length
    """
    bases = ['A', 'C', 'G', 'T']
    result = []
    
    # Convert to base-4
    num = index
    while num > 0:
        result.append(bases[num % 4])
        num //= 4
    
    # Reverse to get correct order
    result.reverse()
    
    # Pad with 'A' to reach target length
    if len(result) < length:
        result = ['A'] * (length - len(result)) + result
    elif len(result) > length:
        raise ValueError(f"Index {index} requires more than {length}bp to encode")
    
    return ''.join(result)


def generate_barcode_mapping(whitelist_file, synthetic_length=25):
    """
    Generate one-to-one mapping from original barcodes to synthetic barcodes.
    
    Args:
        whitelist_file: Path to whitelist file (one barcode per line)
        synthetic_length: Length of synthetic barcodes (default 25, max 31)
    
    Returns:
        dict: {original_barcode: synthetic_barcode}
    """
    if synthetic_length > 31:
        raise ValueError(f"Synthetic barcode length {synthetic_length} exceeds STARsolo limit of 31bp")
    
    mapping = {}
    synthetic_set = set()  # Track for collision detection
    
    with open(whitelist_file) as f:
        for idx, line in enumerate(f):
            original = line.strip()
            if not original:
                continue
            
            # Generate synthetic barcode
            synthetic = encode_index_to_dna(idx, synthetic_length)
            
            # Check for collisions (should never happen with sequential encoding)
            if synthetic in synthetic_set:
                raise RuntimeError(f"Collision detected: synthetic barcode {synthetic} already exists!")
            
            mapping[original] = synthetic
            synthetic_set.add(synthetic)
    
    print(f"Generated {len(mapping)} unique barcode mappings", file=sys.stderr)
    print(f"Synthetic barcode length: {synthetic_length}bp", file=sys.stderr)
    
    return mapping


def write_mapping_table(mapping, output_file):
    """
    Write barcode mapping to TSV file.
    
    Args:
        mapping: dict of original -> synthetic barcodes
        output_file: Output TSV path
    """
    with open(output_file, 'w') as f:
        f.write("original_barcode\tsynthetic_barcode\n")
        for original, synthetic in sorted(mapping.items()):
            f.write(f"{original}\t{synthetic}\n")
    
    print(f"Wrote mapping table: {output_file}", file=sys.stderr)


def write_synthetic_whitelist(mapping, output_file):
    """
    Write whitelist containing only synthetic barcodes.
    
    Args:
        mapping: dict of original -> synthetic barcodes
        output_file: Output whitelist path
    """
    with open(output_file, 'w') as f:
        for synthetic in sorted(set(mapping.values())):
            f.write(f"{synthetic}\n")
    
    print(f"Wrote synthetic whitelist: {output_file}", file=sys.stderr)


def remap_fastq(input_fastq, output_fastq, mapping, barcode_start, barcode_length):
    """
    Rewrite FASTQ file replacing original barcodes with synthetic ones.
    
    Args:
        input_fastq: Input FASTQ path (can be .gz)
        output_fastq: Output FASTQ path (will be .gz)
        mapping: dict of original -> synthetic barcodes
        barcode_start: 0-based start position of barcode in read
        barcode_length: Length of original barcode
    """
    # Open input (handle gzip)
    if str(input_fastq).endswith('.gz'):
        infile = gzip.open(input_fastq, 'rt')
    else:
        infile = open(input_fastq, 'r')
    
    # Open output (always gzip)
    outfile = gzip.open(output_fastq, 'wt')
    
    reads_processed = 0
    reads_remapped = 0
    reads_skipped = 0
    unmapped_barcodes = set()
    
    try:
        while True:
            # Read FASTQ record (4 lines)
            header = infile.readline()
            if not header:
                break
            
            seq = infile.readline().rstrip('\n')
            plus = infile.readline()
            qual = infile.readline().rstrip('\n')
            
            reads_processed += 1
            
            # Extract barcode from sequence
            if len(seq) < barcode_start + barcode_length:
                print(f"Warning: Read {reads_processed} too short ({len(seq)}bp), skipping", file=sys.stderr)
                reads_skipped += 1
                continue
            
            original_barcode = seq[barcode_start:barcode_start + barcode_length]
            
            # Look up synthetic barcode
            if original_barcode not in mapping:
                if original_barcode not in unmapped_barcodes:
                    unmapped_barcodes.add(original_barcode)
                    if len(unmapped_barcodes) <= 10:
                        print(f"Warning: Barcode {original_barcode} not in mapping", file=sys.stderr)
                reads_skipped += 1
                continue
            
            synthetic_barcode = mapping[original_barcode]
            
            # Replace barcode in sequence and quality
            new_seq = seq[:barcode_start] + synthetic_barcode + seq[barcode_start + barcode_length:]
            new_qual = qual[:barcode_start] + qual[barcode_start:barcode_start + len(synthetic_barcode)] + qual[barcode_start + barcode_length:]
            
            # Adjust quality if synthetic is shorter
            if len(synthetic_barcode) < barcode_length:
                # Keep quality scores for synthetic barcode length
                prefix_qual = qual[:barcode_start]
                bc_qual = qual[barcode_start:barcode_start + len(synthetic_barcode)]
                suffix_qual = qual[barcode_start + barcode_length:]
                new_qual = prefix_qual + bc_qual + suffix_qual
            
            # Write remapped record
            outfile.write(header)
            outfile.write(new_seq + '\n')
            outfile.write(plus)
            outfile.write(new_qual + '\n')
            
            reads_remapped += 1
            
            if reads_processed % 100000 == 0:
                print(f"Processed {reads_processed} reads, remapped {reads_remapped}", file=sys.stderr)
    
    finally:
        infile.close()
        outfile.close()
    
    print(f"\nFinal stats:", file=sys.stderr)
    print(f"  Total reads: {reads_processed}", file=sys.stderr)
    print(f"  Remapped: {reads_remapped}", file=sys.stderr)
    print(f"  Skipped: {reads_skipped}", file=sys.stderr)
    print(f"  Unmapped barcodes: {len(unmapped_barcodes)}", file=sys.stderr)
    
    if reads_remapped == 0:
        raise RuntimeError("No reads were successfully remapped!")


def main():
    parser = argparse.ArgumentParser(
        description="Remap 36bp barcodes to synthetic ≤31bp barcodes for STARsolo",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--whitelist', required=True, help='Input whitelist file (original barcodes)')
    parser.add_argument('--fastq', required=True, help='Input FASTQ file with original barcodes')
    parser.add_argument('--output-fastq', required=True, help='Output FASTQ file with synthetic barcodes (.gz)')
    parser.add_argument('--output-whitelist', required=True, help='Output whitelist with synthetic barcodes')
    parser.add_argument('--output-mapping', required=True, help='Output TSV mapping file')
    parser.add_argument('--barcode-start', type=int, default=0, help='0-based barcode start position (default: 0)')
    parser.add_argument('--barcode-length', type=int, default=36, help='Original barcode length (default: 36)')
    parser.add_argument('--synthetic-length', type=int, default=25, help='Synthetic barcode length (default: 25, max: 31)')
    
    args = parser.parse_args()
    
    # Validate
    if args.synthetic_length > 31:
        parser.error(f"Synthetic length {args.synthetic_length} exceeds STARsolo limit of 31bp")
    
    print(f"=== Barcode Remapping for STARsolo ===", file=sys.stderr)
    print(f"Input whitelist: {args.whitelist}", file=sys.stderr)
    print(f"Input FASTQ: {args.fastq}", file=sys.stderr)
    print(f"Original barcode: {args.barcode_length}bp at position {args.barcode_start}", file=sys.stderr)
    print(f"Synthetic barcode: {args.synthetic_length}bp", file=sys.stderr)
    print(f"", file=sys.stderr)
    
    # Step 1: Generate mapping
    print("Step 1: Generating barcode mapping...", file=sys.stderr)
    mapping = generate_barcode_mapping(args.whitelist, args.synthetic_length)
    
    # Step 2: Write mapping table
    print("Step 2: Writing mapping table...", file=sys.stderr)
    write_mapping_table(mapping, args.output_mapping)
    
    # Step 3: Write synthetic whitelist
    print("Step 3: Writing synthetic whitelist...", file=sys.stderr)
    write_synthetic_whitelist(mapping, args.output_whitelist)
    
    # Step 4: Remap FASTQ
    print("Step 4: Remapping FASTQ reads...", file=sys.stderr)
    remap_fastq(args.fastq, args.output_fastq, mapping, args.barcode_start, args.barcode_length)
    
    print("\n=== Remapping complete! ===", file=sys.stderr)


if __name__ == '__main__':
    main()
