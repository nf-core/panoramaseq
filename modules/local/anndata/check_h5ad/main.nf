process ANNDATA_CHECKH5AD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata:0.10.9--d13580e4b297da7c':
        'community.wave.seqera.io/library/anndata:0.10.9--1eab54e300e1e584' }"

    input:
    tuple val(meta), path(h5ad_file)

    output:
    tuple val(meta), path(h5ad_file), emit: h5ad
    path "*.log"                    , emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    check_h5ad.py ${h5ad_file} > ${prefix}_h5ad_validation.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        numpy: \$(python3 -c 'import numpy as np; print(np.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "H5AD validation check completed successfully" > ${prefix}_h5ad_validation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata as ad; print(ad.__version__)')
        numpy: \$(python3 -c 'import numpy as np; print(np.__version__)')
    END_VERSIONS
    """
}
