process STARSOLO_HEATMAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/anndata_numpy_pandas_pillow_pruned:ecb1c72cc96f5cd4' :
        'community.wave.seqera.io/library/anndata_numpy_pandas_pillow_pruned:ecb1c72cc96f5cd4' }"

    input:
    tuple val(meta), path(star_dir)         // STARsolo Solo.out directory
    tuple val(meta2), path(mapping)         // Barcode mapping TSV
    path coords                              // Spatial coordinates CSV

    output:
    tuple val(meta), path("*.png")          , emit: heatmap
    tuple val(meta), path("*_data.tsv")     , emit: data
    tuple val(meta), path("*_stats.json")   , emit: stats
    tuple val(meta), path("*.h5ad")         , optional: true, emit: h5ad
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metric = task.ext.metric ?: 'total_counts'
    def colormap = task.ext.colormap ?: 'viridis'
    def percentile = task.ext.percentile ?: 99.5
    def spot_size = task.ext.spot_size ?: ''
    def save_h5ad = task.ext.save_h5ad ? '--output-h5ad' : ''
    
    """
    # Find the matrix files using symlink-following find
    MATRIX=\$(find -L ${star_dir} -name "matrix.mtx.gz" -path "*/Gene/raw/*" | head -n 1)
    BARCODES=\$(find -L ${star_dir} -name "barcodes.tsv.gz" -path "*/Gene/raw/*" | head -n 1)
    FEATURES=\$(find -L ${star_dir} -name "features.tsv.gz" -path "*/Gene/raw/*" | head -n 1)
    
    if [ -z "\$MATRIX" ] || [ -z "\$BARCODES" ] || [ -z "\$FEATURES" ]; then
        echo "ERROR: Could not find STARsolo output files in ${star_dir}/Gene/raw/"
        exit 1
    fi
    
    # Generate heatmaps for bin sizes: 1 (unbinned), 10, and 50
    for BIN_SIZE in 1 10 50; do
        OUTPUT_BASE="${prefix}_heatmap_bin\${BIN_SIZE}_${metric}_${colormap}"
        
        echo "=== Generating heatmap with bin size \${BIN_SIZE} ==="
        
        python3 ${projectDir}/bin/starsolo_umi_heatmap.py \\
            --matrix \$MATRIX \\
            --barcodes \$BARCODES \\
            --features \$FEATURES \\
            --mapping ${mapping} \\
            --coords ${coords} \\
            --output-png \${OUTPUT_BASE}.png \\
            --output-data \${OUTPUT_BASE}_data.tsv \\
            --output-stats \${OUTPUT_BASE}_stats.json \\
            --sample-id ${prefix} \\
            --bin-size \${BIN_SIZE} \\
            --metric ${metric} \\
            --colormap ${colormap} \\
            --percentile ${percentile} \\
            ${spot_size ? "--spot-size ${spot_size}" : ''} \\
            ${save_h5ad ? "\${save_h5ad} \${OUTPUT_BASE}.h5ad" : ''} \\
            ${args}
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        anndata: \$(python3 -c "import anndata; print(anndata.__version__)")
        pillow: \$(python3 -c "import PIL; print(PIL.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metric = task.ext.metric ?: 'total_counts'
    def colormap = task.ext.colormap ?: 'viridis'
    """
    # Generate stub files for all three bin sizes
    for BIN_SIZE in 1 10 50; do
        OUTPUT_BASE="${prefix}_heatmap_bin\${BIN_SIZE}_${metric}_${colormap}"
        touch \${OUTPUT_BASE}.png
        touch \${OUTPUT_BASE}_data.tsv
        echo '{"sample_id": "${prefix}", "bin_size": '\${BIN_SIZE}'}' > \${OUTPUT_BASE}_stats.json
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        anndata: \$(python3 -c "import anndata; print(anndata.__version__)" 2>/dev/null || echo "0.10.9")
        pillow: \$(python3 -c "import PIL; print(PIL.__version__)" 2>/dev/null || echo "10.0.0")
    END_VERSIONS
    """
}
