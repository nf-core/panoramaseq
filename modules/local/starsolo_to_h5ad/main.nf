process STARSOLO_TO_H5AD {
    tag "${meta.id}"
    label 'process_medium'
    conda "conda-forge::anndata=0.10.9 conda-forge::python=3.11 conda-forge::pandas=2.0.3 conda-forge::scipy=1.11.2 conda-forge::numpy=1.24.3"
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
    # Find STARsolo output files
    MATRIX=\$(find ${star_out_dir} -name "matrix.mtx" -path "*/Solo.out/${feature_type}/filtered/*" | head -1)
    BARCODES=\$(find ${star_out_dir} -name "barcodes.tsv" -path "*/Solo.out/${feature_type}/filtered/*" | head -1)
    FEATURES=\$(find ${star_out_dir} -name "features.tsv" -path "*/Solo.out/${feature_type}/filtered/*" | head -1)
    
    if [ -z "\$MATRIX" ] || [ -z "\$BARCODES" ] || [ -z "\$FEATURES" ]; then
        echo "ERROR: Could not find STARsolo output files in ${star_out_dir}"
        echo "Looking for: Solo.out/${feature_type}/filtered/{matrix.mtx,barcodes.tsv,features.tsv}"
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
