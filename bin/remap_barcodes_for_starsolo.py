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
import csv
import re
import subprocess
import shutil
from pathlib import Path


# IMPROVEMENT 1: Check for pigz availability and use it for parallel compression if available
_PIGZ_AVAILABLE = None

def check_pigz_available():
    """Check if pigz is available on the system PATH."""
    global _PIGZ_AVAILABLE
    if _PIGZ_AVAILABLE is None:
        _PIGZ_AVAILABLE = shutil.which('pigz') is not None
    return _PIGZ_AVAILABLE


def open_fastq_read(filepath):
    """
    Open FASTQ file for reading with pigz if available, otherwise use gzip.
    Returns a file-like object suitable for text reading.
    """
    filepath = str(filepath)
    if filepath.endswith('.gz'):
        if check_pigz_available():
            # Use pigz for parallel decompression
            proc = subprocess.Popen(
                ['pigz', '-dc', filepath],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            return proc.stdout
        else:
            # Fall back to Python gzip
            return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')


def open_fastq_write(filepath):
    """
    Open FASTQ file for writing with pigz if available, otherwise use gzip.
    Returns a file-like object suitable for text writing.
    """
    filepath = str(filepath)
    if filepath.endswith('.gz'):
        if check_pigz_available():
            # Use pigz for parallel compression
            proc = subprocess.Popen(
                ['pigz', '-c'],
                stdin=subprocess.PIPE,
                stdout=open(filepath, 'wb'),
                stderr=subprocess.PIPE,
                text=True
            )
            return proc.stdin
        else:
            # Fall back to Python gzip
            return gzip.open(filepath, 'wt')
    else:
        return open(filepath, 'w')


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
    IMPROVEMENT 6: Removed unnecessary sort - synthetic barcodes are unique by construction.
    
    Args:
        mapping: dict of original -> synthetic barcodes
        output_file: Output whitelist path
    """
    with open(output_file, 'w') as f:
        # Use set to ensure uniqueness (guaranteed by construction but defensive)
        for synthetic in set(mapping.values()):
            f.write(f"{synthetic}\n")
    
    print(f"Wrote synthetic whitelist: {output_file}", file=sys.stderr)


def create_synthetic_coordinates(mapping, coords_file, output_file):
    """
    Create coordinate file with synthetic barcodes instead of original barcodes.
    This allows heatmap generation to skip the reverse mapping step.
    
    Args:
        mapping: dict of original -> synthetic barcodes
        coords_file: Input coordinate CSV file (cell/barcode, x, y)
        output_file: Output coordinate CSV file with synthetic barcodes
    """
    print(f"Creating synthetic barcode coordinates...", file=sys.stderr)
    
    # Read coordinates file using csv module
    with open(coords_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
        # Determine barcode column name
        if not rows:
            print(f"  Warning: Empty coordinate file", file=sys.stderr)
            # Create empty output file
            with open(output_file, 'w') as out_f:
                writer = csv.writer(out_f)
                writer.writerow(['cell', 'x', 'y'])
            return
        
        # Get column names
        fieldnames = reader.fieldnames
        if 'cell' in fieldnames:
            barcode_col = 'cell'
        elif 'barcode' in fieldnames:
            barcode_col = 'barcode'
        else:
            # Assume first column is barcode
            barcode_col = fieldnames[0]
    
    print(f"  Input coordinates: {len(rows)} spots", file=sys.stderr)
    
    # Map original barcodes to synthetic and write output
    n_mapped = 0
    n_unmapped = 0
    
    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(['cell', 'x', 'y'])
        
        for row in rows:
            original_bc = row[barcode_col]
            
            # Map to synthetic barcode
            if original_bc in mapping:
                synthetic_bc = mapping[original_bc]
                writer.writerow([synthetic_bc, row['x'], row['y']])
                n_mapped += 1
            else:
                n_unmapped += 1
    
    print(f"  Mapped {n_mapped}/{len(rows)} spots to synthetic barcodes", file=sys.stderr)
    if n_unmapped > 0:
        print(f"  Warning: {n_unmapped} spots had no matching synthetic barcode", file=sys.stderr)
    print(f"Wrote synthetic coordinates: {output_file}", file=sys.stderr)


def remap_fastq(input_fastq, output_fastq, mapping, barcode_start, barcode_length, umi_length=10, r2_input=None, r2_output=None):
    """
    Rewrite FASTQ file replacing original barcodes with synthetic ones.
    IMPROVEMENT 2: Process R1 and R2 simultaneously to eliminate index set and halve I/O.
    IMPROVEMENT 3: Use zip iterator for efficient FASTQ record reading.
    IMPROVEMENT 4: Assign constant Phred 40 quality (='I') to synthetic barcodes from headers.
    IMPROVEMENT 5: Use regex for robust QUIK header parsing.
    IMPROVEMENT 7: Batch output writes to reduce system call overhead.
    
    Args:
        input_fastq: Input R1 FASTQ path
        output_fastq: Output R1 FASTQ path
        mapping: dict of original -> synthetic barcodes
        barcode_start: 0-based start position of barcode in read
        barcode_length: Length of original barcode
        umi_length: Length of UMI (default: 10bp)
        r2_input: Optional R2 input FASTQ path
        r2_output: Optional R2 output FASTQ path
    """
    # IMPROVEMENT 5: Compile regex for QUIK header parsing once
    # Matches: _calledidx_<digits>_<barcode>$ where barcode is all ACGT
    quik_header_pattern = re.compile(r'_calledidx_\d+_([ACGT]+)$')
    
    # Get synthetic barcode length
    synthetic_length = len(next(iter(mapping.values()))) if mapping else 25
    expected_output_length = synthetic_length + umi_length
    
    # IMPROVEMENT 4: Constant high-quality string for synthetic barcodes
    phred40_quality = 'I' * synthetic_length  # Phred 40 = ASCII 'I'
    
    print(f"Synthetic barcode length: {synthetic_length}bp", file=sys.stderr)
    print(f"UMI length: {umi_length}bp", file=sys.stderr)
    print(f"Output read length will be trimmed to: {expected_output_length}bp", file=sys.stderr)
    if check_pigz_available():
        print(f"Using pigz for parallel compression/decompression", file=sys.stderr)
    
    # IMPROVEMENT 1: Use pigz-aware file openers
    r1_in = open_fastq_read(input_fastq)
    r1_out = open_fastq_write(output_fastq)
    r2_in = open_fastq_read(r2_input) if r2_input else None
    r2_out = open_fastq_write(r2_output) if r2_output else None
    
    # Statistics counters
    reads_processed = 0
    reads_remapped = 0
    reads_skipped = 0
    unmapped_barcodes = set()
    reads_header_extracted = 0
    reads_seq_extracted = 0
    
    # IMPROVEMENT 7: Batch output buffers
    r1_buffer = []
    r2_buffer = []
    batch_size = 10000
    
    try:
        # IMPROVEMENT 3: Use zip iterator for FASTQ record reading
        r1_iter = iter(r1_in)
        r2_iter = iter(r2_in) if r2_in else None
        
        # IMPROVEMENT 2: Process R1 and R2 in lockstep
        for r1_header, r1_seq, r1_plus, r1_qual in zip(r1_iter, r1_iter, r1_iter, r1_iter):
            # Strip newlines
            r1_seq = r1_seq.rstrip('\n')
            r1_qual = r1_qual.rstrip('\n')
            
            # Read corresponding R2 record if present
            r2_record = None
            if r2_iter:
                try:
                    r2_header = next(r2_iter)
                    r2_seq = next(r2_iter)
                    r2_plus = next(r2_iter)
                    r2_qual = next(r2_iter)
                    r2_record = (r2_header, r2_seq, r2_plus, r2_qual)
                except StopIteration:
                    print(f"Warning: R2 ended before R1 at read {reads_processed + 1}", file=sys.stderr)
                    break
            
            reads_processed += 1
            
            # Check minimum read length: need at least umi_length characters to extract UMI.
            # Barcode indels can shorten R1 below (barcode_length + umi_length); the barcode
            # is always taken from the QUIK header, so only the UMI tail is required in seq.
            if len(r1_seq) < barcode_start + umi_length:
                if reads_processed <= 10:
                    print(f"Warning: Read {reads_processed} too short ({len(r1_seq)}bp), skipping", file=sys.stderr)
                reads_skipped += 1
                continue
            
            # IMPROVEMENT 5: Extract barcode using regex (robust to underscores in read names)
            original_barcode = None
            extracted_from_header = False
            
            if r1_header.startswith('@'):
                match = quik_header_pattern.search(r1_header)
                if match:
                    original_barcode = match.group(1)
                    if len(original_barcode) == barcode_length:
                        extracted_from_header = True
            
            # Fall back to sequence extraction if header parsing fails
            if not extracted_from_header:
                original_barcode = r1_seq[barcode_start:barcode_start + barcode_length]
                reads_seq_extracted += 1
            else:
                reads_header_extracted += 1
            
            # Look up synthetic barcode
            if original_barcode not in mapping:
                if original_barcode not in unmapped_barcodes:
                    unmapped_barcodes.add(original_barcode)
                    if len(unmapped_barcodes) <= 10:
                        print(f"Warning: Barcode {original_barcode} not in mapping", file=sys.stderr)
                reads_skipped += 1
                continue
            
            synthetic_barcode = mapping[original_barcode]
            
            # Extract UMI from the tail of the read.
            # When the barcode contains deletions the actual mutated sequence is shorter than
            # barcode_length, so the UMI starts before position barcode_length.  Taking the
            # last umi_length characters is always correct because R1 = mutated_BC + UMI.
            umi = r1_seq[-umi_length:]
            umi_qual = r1_qual[-umi_length:]
            
            # Build new R1 sequence and quality
            new_r1_seq = synthetic_barcode + umi
            
            # IMPROVEMENT 4: Use constant high quality for synthetic barcode if from header
            if extracted_from_header:
                new_r1_qual = phred40_quality + umi_qual
            else:
                # Use original quality if barcode came from sequence
                bc_qual = r1_qual[barcode_start:barcode_start + synthetic_length]
                new_r1_qual = bc_qual + umi_qual
            
            # Verify output length
            if len(new_r1_seq) != expected_output_length:
                if reads_processed <= 10:
                    print(f"Warning: Read {reads_processed} output length mismatch: {len(new_r1_seq)} != {expected_output_length}", file=sys.stderr)
                reads_skipped += 1
                continue
            
            # IMPROVEMENT 7: Append to output buffers
            r1_buffer.append(r1_header)
            r1_buffer.append(new_r1_seq + '\n')
            r1_buffer.append(r1_plus)
            r1_buffer.append(new_r1_qual + '\n')
            
            if r2_record:
                r2_buffer.append(r2_record[0])
                r2_buffer.append(r2_record[1])
                r2_buffer.append(r2_record[2])
                r2_buffer.append(r2_record[3])
            
            reads_remapped += 1
            
            # IMPROVEMENT 7: Flush buffers when batch size reached
            if len(r1_buffer) >= batch_size * 4:
                r1_out.write(''.join(r1_buffer))
                r1_buffer.clear()
                if r2_out:
                    r2_out.write(''.join(r2_buffer))
                    r2_buffer.clear()
            
            # Progress logging
            if reads_processed % 100000 == 0:
                print(f"Processed {reads_processed} reads, remapped {reads_remapped}", file=sys.stderr)
        
        # Flush remaining buffered data
        if r1_buffer:
            r1_out.write(''.join(r1_buffer))
        if r2_buffer and r2_out:
            r2_out.write(''.join(r2_buffer))
    
    finally:
        r1_in.close()
        r1_out.close()
        if r2_in:
            r2_in.close()
        if r2_out:
            r2_out.close()
    
    print(f"\nFinal stats:", file=sys.stderr)
    print(f"  Total reads: {reads_processed}", file=sys.stderr)
    print(f"  Remapped: {reads_remapped}", file=sys.stderr)
    print(f"  Skipped: {reads_skipped}", file=sys.stderr)
    print(f"  Barcodes from QUIK header: {reads_header_extracted}", file=sys.stderr)
    print(f"  Barcodes from sequence: {reads_seq_extracted}", file=sys.stderr)
    print(f"  Unmapped barcodes: {len(unmapped_barcodes)}", file=sys.stderr)
    print(f"  Output read length: {expected_output_length}bp", file=sys.stderr)
    
    if reads_remapped == 0:
        raise RuntimeError("No reads were successfully remapped!")
    
    # IMPROVEMENT 2: R2 is already processed, no need to return passing indices
    if r2_output:
        print(f"  R2 reads written: {reads_remapped}", file=sys.stderr)


# IMPROVEMENT 2: filter_r2_fastq() function removed - R2 filtering now integrated
# into remap_fastq() for simultaneous processing, eliminating the need for
# a separate pass and the memory-intensive passing_indices set.


def main():
    parser = argparse.ArgumentParser(
        description="Remap 36bp barcodes to synthetic ≤31bp barcodes for STARsolo",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--whitelist', required=True, help='Input whitelist file (original barcodes)')
    parser.add_argument('--fastq', required=True, help='Input R1 FASTQ file with original barcodes')
    parser.add_argument('--fastq-r2', required=False, help='Input R2 FASTQ file (optional, will be filtered to match R1)')
    parser.add_argument('--coords', required=False, help='Input barcode coordinates CSV (optional, for generating synthetic coords)')
    parser.add_argument('--output-fastq', required=True, help='Output R1 FASTQ file with synthetic barcodes (.gz)')
    parser.add_argument('--output-fastq-r2', required=False, help='Output R2 FASTQ file (.gz, required if --fastq-r2 provided)')
    parser.add_argument('--output-whitelist', required=True, help='Output whitelist with synthetic barcodes')
    parser.add_argument('--output-mapping', required=True, help='Output TSV mapping file')
    parser.add_argument('--output-coords', required=False, help='Output coordinates CSV with synthetic barcodes (required if --coords provided)')
    parser.add_argument('--barcode-start', type=int, default=0, help='0-based barcode start position (default: 0)')
    parser.add_argument('--barcode-length', type=int, default=36, help='Original barcode length (default: 36)')
    parser.add_argument('--umi-length', type=int, default=10, help='UMI length (default: 10)')
    parser.add_argument('--synthetic-length', type=int, default=25, help='Synthetic barcode length (default: 25, max: 31)')
    
    args = parser.parse_args()
    
    # Validate R2 arguments
    if args.fastq_r2 and not args.output_fastq_r2:
        parser.error("--output-fastq-r2 is required when --fastq-r2 is provided")
    if args.output_fastq_r2 and not args.fastq_r2:
        parser.error("--fastq-r2 is required when --output-fastq-r2 is provided")
    
    # Validate coordinates arguments
    if args.coords and not args.output_coords:
        parser.error("--output-coords is required when --coords is provided")
    if args.output_coords and not args.coords:
        parser.error("--coords is required when --output-coords is provided")
    
    # Validate
    if args.synthetic_length > 31:
        parser.error(f"Synthetic length {args.synthetic_length} exceeds STARsolo limit of 31bp")
    
    print(f"=== Barcode Remapping for STARsolo ===", file=sys.stderr)
    print(f"Input whitelist: {args.whitelist}", file=sys.stderr)
    print(f"Input FASTQ: {args.fastq}", file=sys.stderr)
    print(f"Original barcode: {args.barcode_length}bp at position {args.barcode_start}", file=sys.stderr)
    print(f"UMI length: {args.umi_length}bp", file=sys.stderr)
    print(f"Synthetic barcode: {args.synthetic_length}bp", file=sys.stderr)
    print(f"Output read length: {args.synthetic_length + args.umi_length}bp (trimmed)", file=sys.stderr)
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
    
    # Step 3.5: Create synthetic coordinates if provided
    if args.coords:
        print("Step 3.5: Creating synthetic coordinates...", file=sys.stderr)
        create_synthetic_coordinates(mapping, args.coords, args.output_coords)
    
    # Step 4: Remap FASTQ (R1 and optionally R2 simultaneously)
    if args.fastq_r2:
        print("Step 4: Remapping R1 and R2 FASTQ reads simultaneously...", file=sys.stderr)
    else:
        print("Step 4: Remapping R1 FASTQ reads...", file=sys.stderr)
    
    # IMPROVEMENT 2: Pass R2 files to remap_fastq for simultaneous processing
    remap_fastq(
        args.fastq, 
        args.output_fastq, 
        mapping, 
        args.barcode_start, 
        args.barcode_length, 
        args.umi_length,
        r2_input=args.fastq_r2,
        r2_output=args.output_fastq_r2
    )
    
    print("\n=== Remapping complete! ===", file=sys.stderr)


if __name__ == '__main__':
    main()
