process QUIK_STARSOLO {
    tag "${meta.id}"
    label 'use_gpu'
    // NOTE: Docker is NOT supported. The container is a Singularity .sif image
    // pushed via the ORAS protocol (singularity push oras://...) and cannot be
    // pulled or executed by Docker. Use -profile singularity or -profile apptainer.
    conda null
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/quik-runtime-compile:runtime-compile-v1' :
        error('QUIK_STARSOLO requires Singularity or Apptainer. The container is a native .sif image stored as an OCI artifact (oras://) and cannot be used with Docker. Please run with -profile singularity or -profile apptainer.') }"
    
    input:
    tuple val(meta), path(reads)
    path barcode_file
    
    output:
    tuple val(meta), path("*_R1_filtered.fastq.gz"),                          emit: r1
    tuple val(meta), path("*_R2_filtered.fastq.gz"),                          emit: r2
    tuple val(meta), path("*_whitelist.txt"),                                  emit: whitelist
    tuple val(meta), path("*_barcode_calling_stats.txt"),                      emit: stats
    tuple val(meta), path("*_R1_rejected.fastq.gz"), optional: true,           emit: r1_rejected
    tuple val(meta), path("*_R2_rejected.fastq.gz"), optional: true,           emit: r2_rejected
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Derive barcode/UMI parameters from STARsolo config and read_structure
    def barcode_length = params.starsolo_cb_len
    def umi_length = params.starsolo_umi_len
    def read_structure = params.read_structure
    
    // Calculate barcode start position based on read structure (0-indexed for QUIK)
    //   - BC_UMI: barcode at position 0
    //   - UMI_BC: barcode at position umi_length (after UMI)
    def barcode_start = (read_structure == 'BC_UMI') ? 0 : umi_length
    
    // QUIK-specific parameters
    def strategy = params.strategy
    def distance_measure = params.distance_measure
    def rejection_threshold = params.rejection_threshold ?: (barcode_length * 0.25).toInteger()
    
    """
    # Using runtime-compiled QUIK for flexible barcode lengths
    echo "=== QUIK Runtime Compilation for STARsolo Workflow ==="
    echo "Building QUIK executable..."
    echo "Parameters: barcode_length=${barcode_length}, rejection_threshold=${rejection_threshold}"
    echo "Strategy: ${strategy}, Distance: ${distance_measure}, Barcode start: ${barcode_start}"
    echo "====================================="
    
    # Step 1: Compile QUIK (in writable /tmp directory)
    mkdir -p /tmp/quik_build_${prefix} && cd /tmp/quik_build_${prefix}
    cmake /opt/quik
    make -j${task.cpus}
    QUIK_EXEC=\$(pwd)/single_strategy_benchmark_fastq_paired
    cd -
    
    echo "QUIK compiled successfully: \${QUIK_EXEC}"
    
    # Step 2: Decompress input FASTQ files
    echo "Decompressing input FASTQ files..."
    gunzip -c ${reads[0]} > input_R1.fastq
    gunzip -c ${reads[1]} > input_R2.fastq
    
    # Step 3: Extract just the barcode sequences from the CSV file
    echo "Extracting barcode sequences..."
    tail -n +2 ${barcode_file} | cut -d',' -f1 > barcodes_only.txt
    
    # Step 4: Run QUIK barcode calling with flexible parameters
    \${QUIK_EXEC} \\
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
    
    # Step 5: Use reference barcodes as whitelist (not observed sequences)
    # The whitelist should contain the reference barcodes that QUIK matched against,
    # not the observed sequences with sequencing errors
    echo "Creating whitelist from reference barcodes..."
    cp barcodes_only.txt ${prefix}_whitelist.txt
    
    # Step 6: Validate whitelist
    whitelist_count=\$(wc -l < ${prefix}_whitelist.txt)
    echo "Whitelist contains \${whitelist_count} reference barcodes"
    
    # Step 7: Compress output FASTQ files
    gzip ${prefix}_R1_filtered.fastq
    gzip ${prefix}_R2_filtered.fastq

    # Step 8: Extract rejected reads (optional, for Columba barcode rescue)
    if [ "${params.enable_columba_rescue}" = "true" ]; then
        echo "Extracting rejected reads for Columba barcode rescue..."
        
        # Extract read IDs from filtered FASTQ files (every 4th line starting from line 1)
        # QUIK appends _calledidx_<idx>_<barcode> to read IDs, so we extract only the original read ID
        echo "Extracting filtered read IDs..."
        zcat ${prefix}_R1_filtered.fastq.gz | awk 'NR % 4 == 1 {
            read_id = substr(\$1, 2)
            split(read_id, parts, "_calledidx_")
            print parts[1]
        }' > filtered_ids.txt
        
        filtered_count=\$(wc -l < filtered_ids.txt)
        echo "Filtered read count: \${filtered_count}"
        
        # Extract rejected reads: find reads NOT in the filtered IDs list
        # Using FNR % 4 logic to properly parse FASTQ records without duplicates
        echo "Extracting rejected R1 reads..."
        awk 'NR==FNR{ids[\$1]=1; next} 
             FNR % 4 == 1 {
                 hdr = \$0
                 read_id = substr(\$1, 2)
                 split(read_id, parts, " ")
                 id = parts[1]
             }
             FNR % 4 == 2 {seq = \$0}
             FNR % 4 == 3 {plus = \$0}
             FNR % 4 == 0 {
                 qual = \$0
                 if (!(id in ids)) {
                     print hdr"\\n"seq"\\n"plus"\\n"qual
                 }
             }' filtered_ids.txt input_R1.fastq | gzip > ${prefix}_R1_rejected.fastq.gz
        
        echo "Extracting rejected R2 reads..."
        awk 'NR==FNR{ids[\$1]=1; next} 
             FNR % 4 == 1 {
                 hdr = \$0
                 read_id = substr(\$1, 2)
                 split(read_id, parts, " ")
                 id = parts[1]
             }
             FNR % 4 == 2 {seq = \$0}
             FNR % 4 == 3 {plus = \$0}
             FNR % 4 == 0 {
                 qual = \$0
                 if (!(id in ids)) {
                     print hdr"\\n"seq"\\n"plus"\\n"qual
                 }
             }' filtered_ids.txt input_R2.fastq | gzip > ${prefix}_R2_rejected.fastq.gz
        
        rejected_count=\$(zcat ${prefix}_R1_rejected.fastq.gz | wc -l)
        rejected_reads=\$((rejected_count / 4))
        echo "Rejected R1 reads: \${rejected_reads}"
        
        # Clean up
        rm filtered_ids.txt
    fi

    # Step 9: Clean up
    rm input_R1.fastq input_R2.fastq barcodes_only.txt
    rm -rf /tmp/quik_build_${prefix}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-runtime-flex")
        cuda: \$(nvcc --version 2>/dev/null | grep release | cut -d' ' -f6 | cut -d',' -f1 || echo "12.6.0")
        barcode_length: \$(echo "${barcode_length}")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_length = params.starsolo_cb_len
    """
    touch ${prefix}_R1_filtered.fastq.gz
    touch ${prefix}_R2_filtered.fastq.gz
    touch ${prefix}_whitelist.txt
    touch ${prefix}_barcode_calling_stats.txt
    touch ${prefix}_R1_rejected.fastq.gz
    touch ${prefix}_R2_rejected.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quik: \$(echo "2.0-runtime-flex")
        barcode_length: \$(echo "${barcode_length}")
    END_VERSIONS
    """
}
