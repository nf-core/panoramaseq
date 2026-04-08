process EXTRACT_WHITELIST_FROM_COORDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(coords)  // CSV with columns: cell/barcode, x, y

    output:
    tuple val(meta), path("*_whitelist.txt"), emit: whitelist
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract first column (barcodes) from CSV, skip header
    tail -n +2 ${coords} | cut -d',' -f1 > ${prefix}_whitelist.txt
    
    echo "Extracted \$(wc -l < ${prefix}_whitelist.txt) barcodes to whitelist from ${coords}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_whitelist.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
