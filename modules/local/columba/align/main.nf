process COLUMBA_ALIGN {
    tag "$meta.id"
    label 'process_high'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/columba_vanilla:latest' :
        'quay.io/francoaps/columba_vanilla:latest' }"
    
    input:
    tuple val(meta), path(reads), path(index_files)
    path binaries_dir   // build_Vanilla/ directory from COLUMBA_BUILD
    
    output:
    tuple val(meta), path("*_alignment.sam"), emit: sam
    tuple val(meta), path("*_columba_stats.txt"), emit: stats
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def identity_threshold = params.columba_identity_threshold ?: 83
    def barcode_length = params.starsolo_cb_len ?: 36
    def barcode_start = (params.read_structure == 'UMI_BC') ? (params.starsolo_umi_len ?: 10) : 0
    def barcode_end = barcode_start + barcode_length - 1
    def barcode_window = "${barcode_start}-${barcode_end}"
    def threads = task.cpus

    // Use R1 for rejected reads (contains barcode region)
    def input_fastq = reads instanceof List ? reads[0] : reads
    """
    echo "=== Running Columba Alignment ==="
    echo "Identity threshold: ${identity_threshold}"
    echo "Barcode window: ${barcode_window}  (barcode length: ${barcode_length})"
    echo "Threads: ${threads}"
    echo "=================================="

    # Reference columba binary directly from staged binaries directory
    COLUMBA_BIN="${binaries_dir}/columba"
    chmod +x "\${COLUMBA_BIN}"
    echo "Using columba: \${COLUMBA_BIN}"

    # Run Columba alignment
    START_TIME=\$(date +%s.%N)
    "\${COLUMBA_BIN}" \\
        -f ${input_fastq} \\
        -r ${meta.id}_index \\
        -I ${identity_threshold} \\
        -T ${barcode_window} \\
        -t ${threads} \\
        ${args} \\
        2>&1 | tee ${prefix}_columba.log

    END_TIME=\$(date +%s.%N)
    DURATION=\$(awk -v start=\$START_TIME -v end=\$END_TIME 'BEGIN {printf "%.3f", end-start}')

    # Rename output SAM file
    mv ColumbaOutput.sam ${prefix}_alignment.sam

    # Create statistics file
    {
        echo "Columba Barcode Rescue Statistics"
        echo "=========================================="
        echo "Sample ID: ${meta.id}"
        echo "Identity threshold (-I): ${identity_threshold}"
        echo "Barcode window (-T): ${barcode_window}"
        echo "Threads: ${threads}"
        echo ""
        TOTAL_READS=\$(zcat ${input_fastq} | wc -l | awk '{print \$1/4}')
        ALIGNED_READS=\$(grep -cv '^@' ${prefix}_alignment.sam || echo 0)
        UNMAPPED=\$(awk '\$2 == 4 {count++} END {print count+0}' ${prefix}_alignment.sam)
        MAPPED=\$((ALIGNED_READS - UNMAPPED))
        echo "Alignment statistics:"
        echo "  Total reads: \$TOTAL_READS"
        echo "  Mapped reads: \$MAPPED"
        echo "  Unmapped reads: \$UNMAPPED"
        if [ "\$TOTAL_READS" -gt 0 ]; then
            awk -v mapped=\$MAPPED -v total=\$TOTAL_READS 'BEGIN {printf "  Mapping rate: %.2f%%\\n", (mapped/total)*100}'
        fi
        echo ""
        echo "Timing information:"
        echo "  Total alignment time: \${DURATION} seconds"
        if [ "\$TOTAL_READS" -gt 0 ]; then
            awk -v dur=\$DURATION -v reads=\$TOTAL_READS 'BEGIN {printf "  Time per read: %.6f ms\\n", (dur/reads)*1000}'
        fi
    } > ${prefix}_columba_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba: "\$("\${COLUMBA_BIN}" --version 2>&1 | head -1 || echo 'unknown')"
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_alignment.sam
    touch ${prefix}_columba_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba: "stub-version"
    END_VERSIONS
    """
}
