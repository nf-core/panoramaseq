process ANNDATA_MAKEH5AD_SINGLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata:0.10.9--d13580e4b297da7c':
        'community.wave.seqera.io/library/anndata:0.10.9--1eab54e300e1e584' }"

    input:
    tuple val(meta), path(count_tsv), path(coords_csv)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name_arg = meta.id ? "--sample-name ${meta.id}" : ""
    """
    # Check if input TSV file is empty (only header or completely empty)
    echo "Checking input file ${count_tsv} for data..."
    lines=\$(zcat ${count_tsv} | wc -l)
    echo "File has \$lines lines (including header)"

    # If file has 1 or fewer lines (header only or empty), exit with error
    if [ "\$lines" -le 1 ]; then
        echo ""
        echo "ERROR: Input count file ${count_tsv} is empty (only \$lines lines)"
        echo "Sample: ${meta.id}"
        echo ""
        echo "This typically happens when:"
        echo "  1. The GTF file has no features overlapping with aligned reads for this sample"
        echo "  2. Test data is too small (e.g., test profile with truncated GTF)"
        echo "  3. Barcode calling filtered out all reads for this sample"
        echo "  4. UMI counting produced no results for this sample"
        echo ""
        echo "NOTE: This sample will be skipped. Other samples may still succeed."
        echo ""
        echo "Suggestions:"
        echo "  - Use a complete GTF file matching your reference genome"
        echo "  - Check UMICOUNT logs in the work directory for warnings"
        echo "  - Verify barcode calling produced valid output files"
        echo "  - Ensure your test data has reads that map to genes in the GTF"

        # Create a dummy output to satisfy Nextflow output requirements
        # This allows the errorStrategy = 'ignore' to work properly
        touch ${prefix}.h5ad.failed

        exit 1
    fi

    tsv_to_h5ad_single.py \\
        ${count_tsv} \\
        --coords ${coords_csv} \\
        --output ${prefix}.h5ad \\
        ${sample_name_arg} \\
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
