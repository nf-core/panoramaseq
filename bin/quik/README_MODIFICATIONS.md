# QUIK Barcode Calling - Modified Version

## Overview
This is a modified version of the QUIK barcode calling software, adapted for use in Nextflow pipelines with flexible command-line parameters.

## Key Modifications

### 1. New Executable: `single_strategy_benchmark_fastq_paired`

**Purpose**: Process paired-end FASTQ files with configurable barcode calling parameters.

**Command-line interface**:
```bash
single_strategy_benchmark_fastq_paired \
    <barcode_file> \
    <r1_fastq> \
    <r2_fastq> \
    <barcode_start> \
    <barcode_length> \
    <strategy> \
    <distance_measure> \
    <rejection_threshold> \
    <output_r1> \
    <output_r2>
```

**Parameters**:
- `barcode_file`: Text file with one barcode per line
- `r1_fastq`: Forward reads FASTQ file
- `r2_fastq`: Reverse reads FASTQ file
- `barcode_start`: Starting position of barcode in R1 (0-indexed)
- `barcode_length`: Length of barcode sequence
- `strategy`: Calling strategy (e.g., 4_7_mer_gpu_v1, 5_mer_gpu_v4, two_step_gpu_v1)
- `distance_measure`: LEVENSHTEIN, SEQUENCE_LEVENSHTEIN, or PSEUDO_DISTANCE
- `rejection_threshold`: Maximum distance for assignment
- `output_r1`: Output file for filtered R1 reads
- `output_r2`: Output file for filtered R2 reads

### 2. Supported Strategies

**GPU Strategies**:
- `4_7_mer_gpu_v1`, `5_7_mer_gpu_v1`, `6_7_mer_gpu_v1`, `7_7_mer_gpu_v1` (v1 variants)
- `4_mer_gpu_v2`, `5_mer_gpu_v2` (v2 variants)
- `4_mer_gpu_v3`, `5_mer_gpu_v3` (v3 variants)
- `4_mer_gpu_v4`, `5_mer_gpu_v4`, `6_mer_gpu_v4`, `7_mer_gpu_v4` (v4 variants)
- `4_mer_gpu_v5`, `5_mer_gpu_v5` (v5 variants)
- `two_step_gpu_v1` (two-stage calling: 5-mer then 7-mer)

**CPU Strategies**:
- `host_v1`, `host_v2` (5-mer filtering on CPU)
- `two_step_host_v1` (two-stage calling on CPU)

### 3. Critical Bug Fixes

#### Sequence Length Handling
**Problem**: Original code had `assert(seq.size() == LENGTH)` in sequence constructor, causing crashes on variable-length sequences.

**Solution**: 
- Normalize all sequences to `SEQUENCE_LENGTH` by padding with 'A' or truncating
- Extract barcodes from reads with padding for short sequences
- Normalize barcode file sequences to handle length mismatches

#### FASTQ Format Support
**Problem**: Original code only read plain text sequence files.

**Solution**: 
- Implemented full FASTQ parser for paired-end reads
- Reads header, sequence, plus line, and quality scores
- Processes both R1 and R2 files simultaneously

#### Barcode Extraction
**Problem**: No support for extracting barcodes from specific positions in reads.

**Solution**:
- Configurable barcode start position and length
- Automatic padding with 'A' for reads shorter than expected
- Handles edge cases gracefully

### 4. Output Files

The program generates:
1. **Filtered R1 FASTQ**: Contains only reads that were successfully assigned to barcodes
2. **Filtered R2 FASTQ**: Corresponding mate pairs for assigned reads
3. **Statistics (stdout/stderr)**: Detailed statistics including:
   - Total reads processed
   - Assignment vs rejection counts
   - Processing time and speed
   - Barcode distribution (top 10)
   - Short read padding statistics

### 5. Building

**Requirements**:
- CMake 3.22+
- CUDA 12.0+ with compute capability 6.0+ (sm_60, sm_80)
- OpenMP support
- C++14 compiler

**Build commands**:
```bash
mkdir -p build
cd build
cmake ..
make single_strategy_benchmark_fastq_paired
```

### 6. Integration with Nextflow

The modified version is designed to work seamlessly with Nextflow:
- Source code stored in `bin/quik/` directory
- Built on-the-fly in each workflow execution
- No pre-compilation needed
- Uses HPC modules (CUDA/12.6.0, CMake/3.29.3-GCCcore-13.3.0)

## Modified Files

1. **src/single_strategy_benchmark_fastq_paired.cu** (NEW)
   - Main implementation with CLI argument parsing
   - FASTQ file reading and barcode extraction
   - Dynamic strategy selection
   - Sequence normalization and padding

2. **CMakeLists.txt** (MODIFIED)
   - Added new executable target: `single_strategy_benchmark_fastq_paired`
   - Configured with CUDA separable compilation and OpenMP

## Example Usage

```bash
# Load required modules
module load CUDA/12.6.0 CMake/3.29.3-GCCcore-13.3.0

# Build
cd build
cmake ..
make single_strategy_benchmark_fastq_paired

# Run barcode calling
./single_strategy_benchmark_fastq_paired \
    barcodes.txt \
    sample_R1.fastq \
    sample_R2.fastq \
    9 \
    36 \
    4_7_mer_gpu_v1 \
    SEQUENCE_LEVENSHTEIN \
    8 \
    output_R1_filtered.fastq \
    output_R2_filtered.fastq \
    > stats.txt 2>&1
```

## Performance Notes

- GPU acceleration provides significant speedup (typically 0.05-0.5 ms per read)
- Two-step strategies useful for recovering initially unassigned reads
- Higher rejection thresholds increase recall but may reduce precision
- `SEQUENCE_LEVENSHTEIN` distance typically more accurate than `LEVENSHTEIN`

## Version History

**v1.0.0** (October 2025)
- Initial modified version for Nextflow integration
- Added flexible CLI interface
- Implemented FASTQ support
- Fixed sequence length assertion bug
- Added barcode extraction with padding

## Original Source

Based on the QUIK barcode calling algorithm.
See original repository documentation for algorithm details.
