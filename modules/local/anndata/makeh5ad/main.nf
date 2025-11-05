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
    # Check if input TSV files are empty (only header or completely empty)
    echo "Checking input files for data..."
    valid_files=""
    empty_files=""
    empty_count=0
    valid_count=0
    
    for tsv in ${tsv_files}; do
        # Get uncompressed size (number of lines)
        lines=\$(zcat "\$tsv" | wc -l)
        # If file has 1 or fewer lines (header only or empty), mark as empty
        if [ "\$lines" -le 1 ]; then
            echo "WARNING: File \$tsv appears to be empty (only \$lines lines) - skipping"
            empty_files="\${empty_files} \$tsv"
            empty_count=\$((empty_count + 1))
        else
            echo "File \$tsv has \$lines lines (including header) - valid"
            valid_files="\${valid_files} \$tsv"
            valid_count=\$((valid_count + 1))
        fi
    done
    
    echo ""
    echo "Summary: \$valid_count valid files, \$empty_count empty files"
    
    # If ALL files are empty, exit with error
    if [ \$valid_count -eq 0 ]; then
        echo ""
        echo "ERROR: All input count files are empty. Cannot create H5AD file."
        echo "Empty files:"
        echo "\$empty_files"
        echo ""
        echo "This typically happens when:"
        echo "  1. The GTF file has no features overlapping with aligned reads"
        echo "  2. Test data is too small (e.g., test profile with truncated GTF)"
        echo "  3. Barcode calling filtered out all reads"
        echo "  4. UMI counting produced no results"
        echo ""
        echo "Suggestions:"
        echo "  - Use a complete GTF file matching your reference genome"
        echo "  - Check UMICOUNT logs in the work directory for warnings"
        echo "  - Verify barcode calling produced valid output files"
        echo "  - Ensure your test data has reads that map to genes in the GTF"
        exit 1
    fi
    
    # If some files are empty but others are valid, proceed with valid files only
    if [ \$empty_count -gt 0 ]; then
        echo ""
        echo "NOTE: Proceeding with \$valid_count valid files, skipping \$empty_count empty files"
        echo "Empty files skipped: \$empty_files"
        echo ""
    fi
    
    # Process only the valid files
    tsv_to_h5ad.py \\
        \$valid_files \\
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