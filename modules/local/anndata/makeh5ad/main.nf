process ANNDATA_MAKEH5AD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata:0.10.9--d13580e4b297da7c':
        'community.wave.seqera.io/library/anndata:0.10.9--1eab54e300e1e584' }"

    input:
    tuple val(meta), path(count_tsvs), path(coords_csv)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tsv_files = count_tsvs instanceof java.util.Collection ? count_tsvs.join(' ') : count_tsvs
    """
    tsv_to_h5ad.py \\
        ${tsv_files} \\
        --coords ${coords_csv} \\
        --output ${prefix}.h5ad \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
    END_VERSIONS
    """
}