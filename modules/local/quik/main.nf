process QUIK_BARCODE_CALLING {
    tag "${meta.id}"
    label 'gpu_process'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/quik-cuda:latest' :
        'quay.io/francoaps/quik-cuda:latest' }"
    
    input:
    tuple val(meta), path(reads)
    path barcode_file
    
    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*_barcode_calling_stats.txt"), emit: stats
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Extract parameters from params (pipeline-level configuration)
    // These can be overridden via command line: --barcode_start 10 --barcode_length 40 etc.
    def barcode_start = params.barcode_start
    def barcode_length = params.barcode_length
    def strategy = params.strategy
    def distance_measure = params.distance_measure
    def rejection_threshold = params.rejection_threshold
    
    """
    # Debug: Check environment and available tools
    echo "=== Environment Debug ==="
    echo "PATH: \$PATH"
    echo "PWD: \$(pwd)"
    which cmake || echo "cmake not found in PATH"
    which make || echo "make not found in PATH"
    which g++ || echo "g++ not found in PATH"
    which nvcc || echo "nvcc not found in PATH"
    cmake --version || echo "cmake version check failed"
    nvcc --version || echo "nvcc version check failed"
    echo "========================="
    
    # Copy quik source from bin directory
    cp -r ${projectDir}/bin/quik .
    
    # Build the executable using HPC modules
    cd quik
    # Clean any previous build artifacts to avoid cache conflicts
    rm -rf build
    mkdir -p build
    cd build
    
    # Pass pipeline parameters to CMake as compile-time definitions
    echo "Configuring QUIK with SEQUENCE_LENGTH=${barcode_length} and REJECTION_THRESHOLD=${rejection_threshold}"
    cmake -DSEQUENCE_LENGTH=${barcode_length} -DREJECTION_THRESHOLD=${rejection_threshold} ..
    make -j${task.cpus}
    
    # Copy executable to working directory
    echo "Files in build directory:"
    ls -la
    echo "Current directory: \$(pwd)"
    WORK_DIR=\$(pwd | sed 's|/quik/build||')
    echo "Work directory: \$WORK_DIR"
    echo "Copying executable to work directory..."
    cp single_strategy_benchmark_fastq_paired \$WORK_DIR/
    cd ../..
    echo "Files in work directory after copy:"
    ls -la
    echo "Making executable..."
    chmod +x single_strategy_benchmark_fastq_paired
    
    # Decompress input FASTQ files (QUIK requires uncompressed input)
    echo "Decompressing input FASTQ files..."
    gunzip -c ${reads[0]} > input_R1.fastq
    gunzip -c ${reads[1]} > input_R2.fastq
    
    # Extract just the barcode sequences from the CSV file (skip header, take first column)
    echo "Extracting barcode sequences..."
    tail -n +2 ${barcode_file} | cut -d',' -f1 > barcodes_only.txt
    
    # Run quik_clean barcode calling using the built executable
    ./single_strategy_benchmark_fastq_paired \\
        barcodes_only.txt \\
        input_R1.fastq \\
        input_R2.fastq \\
        ${barcode_start} \\
        ${barcode_length} \\
        ${strategy} \\
        ${distance_measure} \\
        ${rejection_threshold} \\
        ${prefix}_R1_filtered.fastq \\
        ${prefix}_R2_filtered.fastq \\
        > ${prefix}_barcode_calling_stats.txt 2>&1
    
    # Clean up decompressed input files to save space
    rm input_R1.fastq input_R2.fastq barcodes_only.txt
    
    # Compress output FASTQ files to match pipeline expectations
    gzip ${prefix}_R1_filtered.fastq
    gzip ${prefix}_R2_filtered.fastq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik_clean: \$(echo "1.0.0")
        cuda: \$(nvcc --version | grep release | cut -d' ' -f6 | cut -d',' -f1)
        cmake: \$(cmake --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_filtered.fastq.gz
    touch ${prefix}_R2_filtered.fastq.gz
    touch ${prefix}_barcode_calling_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik_clean: \$(echo "1.0.0")
    END_VERSIONS
    """
}
