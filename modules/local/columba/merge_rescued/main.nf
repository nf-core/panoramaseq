process MERGE_RESCUED_READS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    // QUIK-filtered reads (keyed by meta, joined before calling this process)
    tuple val(meta), path(r1_filtered),    path(r2_filtered),    path(whitelist_quik),
                     path(r1_rescued),     path(r2_rescued),     path(whitelist_rescued)

    output:
    tuple val(meta), path("*_R1_merged.fastq.gz"),  emit: r1
    tuple val(meta), path("*_R2_merged.fastq.gz"),  emit: r2
    tuple val(meta), path("*_whitelist_merged.txt"), emit: whitelist
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Concatenate gzip streams (valid: gzip supports multi-stream concatenation)
    cat ${r1_filtered} ${r1_rescued} > ${prefix}_R1_merged.fastq.gz
    cat ${r2_filtered} ${r2_rescued} > ${prefix}_R2_merged.fastq.gz

    # Merge whitelists: union of both, deduplicated and sorted
    cat ${whitelist_quik} ${whitelist_rescued} | sort -u > ${prefix}_whitelist_merged.txt

    QUIK_N=\$(wc -l < ${whitelist_quik})
    RESCUED_N=\$(wc -l < ${whitelist_rescued})
    MERGED_N=\$(wc -l < ${prefix}_whitelist_merged.txt)
    echo "QUIK whitelist: \$QUIK_N barcodes"
    echo "Columba rescued: \$RESCUED_N new barcodes"
    echo "Merged whitelist: \$MERGED_N unique barcodes"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_merged.fastq.gz
    touch ${prefix}_R2_merged.fastq.gz
    touch ${prefix}_whitelist_merged.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub-version"
    END_VERSIONS
    """
}
