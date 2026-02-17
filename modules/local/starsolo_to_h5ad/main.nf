process STARSOLO_TO_H5AD {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata:0.10.9--d13580e4b297da7c' :
        'community.wave.seqera.io/library/anndata:0.10.9--1eab54e300e1e584' }"
    
    input:
    tuple val(meta), path(star_out_dir)
    path barcode_coords
    
    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def feature_type = task.ext.feature_type ?: "Gene"
    """
    # Find STARsolo output files (can be gzipped or uncompressed)
    MATRIX=\$(find ${star_out_dir} -path "*/${feature_type}/raw/*" \\( -name "matrix.mtx.gz" -o -name "matrix.mtx" \\) | head -1)
    BARCODES=\$(find ${star_out_dir} -path "*/${feature_type}/raw/*" \\( -name "barcodes.tsv.gz" -o -name "barcodes.tsv" \\) | head -1)
    FEATURES=\$(find ${star_out_dir} -path "*/${feature_type}/raw/*" \\( -name "features.tsv.gz" -o -name "features.tsv" \\) | head -1)
    
    if [ -z "\$MATRIX" ] || [ -z "\$BARCODES" ] || [ -z "\$FEATURES" ]; then
        echo "ERROR: Could not find STARsolo output files in ${star_out_dir}"
        echo "Looking for: ${feature_type}/raw/{matrix.mtx[.gz],barcodes.tsv[.gz],features.tsv[.gz]}"
        exit 1
    fi
    
    echo "Found STARsolo output files:"
    echo "Matrix: \$MATRIX"
    echo "Barcodes: \$BARCODES"
    echo "Features: \$FEATURES"
    
    # Convert to H5AD
    starsolo_to_h5ad.py \\
        --matrix "\$MATRIX" \\
        --barcodes "\$BARCODES" \\
        --features "\$FEATURES" \\
        --coords "${barcode_coords}" \\
        --output "${prefix}.h5ad" \\
        --sample-id "${meta.id}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        anndata: \$(python3 -c "import anndata; print(anndata.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        anndata: 0.10.9
        pandas: 2.0.3
        scipy: 1.11.2
    END_VERSIONS
    """
}
