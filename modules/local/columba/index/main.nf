process COLUMBA_INDEX {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/columba_vanilla:latest' :
        'quay.io/francoaps/columba_vanilla:latest' }"
    
    input:
    tuple val(meta), path(barcode_fasta)
    path binaries_dir   // build_Vanilla/ directory from COLUMBA_BUILD
    
    output:
    tuple val(meta), path("*_index*"), emit: index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Reference columba_build binary directly from staged binaries directory
    COLUMBA_BUILD_BIN="${binaries_dir}/columba_build"
    chmod +x "\${COLUMBA_BUILD_BIN}"
    echo "Using columba_build: \${COLUMBA_BUILD_BIN}"

    # Build Columba index from FASTA file
    "\${COLUMBA_BUILD_BIN}" \\
        -r ${prefix}_index \\
        -f ${barcode_fasta} \\
        ${args}

    # Verify index files were created
    ls -lh ${prefix}_index*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba_build: "\$("\${COLUMBA_BUILD_BIN}" --version 2>&1 | head -1 || echo 'unknown')"
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_index.columba
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba_build: "stub-version"
    END_VERSIONS
    """
}
