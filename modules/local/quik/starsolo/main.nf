process QUIK_STARSOLO {
    tag "${meta.id}"
    label 'use_gpu'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/quik-cuda:prebuilt-36bp-v2' :
        'quay.io/francoaps/quik-cuda:prebuilt-36bp-v2' }"
    
    input:
    tuple val(meta), path(reads)
    path barcode_file
    
    output:
    tuple val(meta), path("*_R1_filtered.fastq.gz"), emit: r1
    tuple val(meta), path("*_R2_filtered.fastq.gz"), emit: r2
    path "*_whitelist.txt", emit: whitelist
    tuple val(meta), path("*_barcode_calling_stats.txt"), emit: stats
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Extract parameters from params
    def barcode_start = params.barcode_start
    def barcode_length = 36  // FIXED in pre-built binary
    def strategy = params.strategy
    def distance_measure = params.distance_measure
    def rejection_threshold = 8  // FIXED in pre-built binary
    
    """
    # Using pre-built QUIK binary for STARsolo workflow
    echo "=== QUIK for STARsolo Workflow ==="
    echo "Binary location: \$(which quik)"
    echo "Configured for: SEQUENCE_LENGTH=36, REJECTION_THRESHOLD=8"
    echo "Runtime parameters: barcode_start=${barcode_start}, strategy=${strategy}, distance=${distance_measure}"
    echo "Output: R1 with corrected barcodes + whitelist for STARsolo"
    echo "===================================="
    
    # Decompress input FASTQ files
    echo "Decompressing input FASTQ files..."
    gunzip -c ${reads[0]} > input_R1.fastq
    gunzip -c ${reads[1]} > input_R2.fastq
    
    # Extract just the barcode sequences from the CSV file
    echo "Extracting barcode sequences..."
    tail -n +2 ${barcode_file} | cut -d',' -f1 > barcodes_only.txt
    
    # Run QUIK barcode calling
    quik \\
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
    
    # Extract unique called barcodes from R1 to create whitelist
    # Note: Barcodes are at positions 11-46 (after 10bp UMI)
    echo "Generating whitelist from called barcodes..."
    awk 'NR%4==2 {print substr(\$0, 11, 36)}' ${prefix}_R1_filtered.fastq | \\
        sort -u > ${prefix}_whitelist.txt
    
    # Validate whitelist
    whitelist_count=\$(wc -l < ${prefix}_whitelist.txt)
    echo "Whitelist contains \${whitelist_count} unique barcodes"
    
    # Clean up decompressed input files
    rm input_R1.fastq input_R2.fastq barcodes_only.txt
    
    # Compress output FASTQ files
    gzip ${prefix}_R1_filtered.fastq
    gzip ${prefix}_R2_filtered.fastq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-prebuilt-36bp")
        cuda: \$(nvcc --version 2>/dev/null | grep release | cut -d' ' -f6 | cut -d',' -f1 || echo "12.6.0")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_filtered.fastq.gz
    touch ${prefix}_R2_filtered.fastq.gz
    touch ${prefix}_whitelist.txt
    touch ${prefix}_barcode_calling_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-prebuilt-36bp")
    END_VERSIONS
    """
}
