process REMAP_BARCODES_FOR_STARSOLO {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/pigz_python:b7952300f0cafe48' :
        'community.wave.seqera.io/library/pigz_python:b7952300f0cafe48' }"
    
    input:
    tuple val(meta), path(reads)        // R1 reads
    tuple val(meta2), path(whitelist)   // Whitelist
    tuple val(meta3), path(reads_r2)    // R2 reads
    path coords                          // Barcode coordinates (optional)
    
    output:
    tuple val(meta), path("*_synthetic.fastq.gz"), emit: reads
    tuple val(meta), path("*_R2_remapped.fastq.gz"), emit: reads_r2
    path "*_whitelist_synthetic.txt"              , emit: whitelist
    tuple val(meta), path("*_barcode_mapping.tsv"), emit: mapping
    path "*_coords_synthetic.csv"                 , optional: true, emit: coords
    path "versions.yml"                           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_start = task.ext.barcode_start ?: 0
    def barcode_length = task.ext.barcode_length ?: 36
    def umi_length = task.ext.umi_length ?: 10
    def synthetic_length = task.ext.synthetic_length ?: 25
    def coords_arg = coords.name != 'NO_FILE' ? "--coords ${coords} --output-coords ${prefix}_coords_synthetic.csv" : ""
    """
    remap_barcodes_for_starsolo.py \\
        --whitelist ${whitelist} \\
        --fastq ${reads} \\
        --fastq-r2 ${reads_r2} \\
        --output-fastq ${prefix}_synthetic.fastq.gz \\
        --output-fastq-r2 ${prefix}_R2_remapped.fastq.gz \\
        --output-whitelist ${prefix}_whitelist_synthetic.txt \\
        --output-mapping ${prefix}_barcode_mapping.tsv \\
        --barcode-start ${barcode_start} \\
        --barcode-length ${barcode_length} \\
        --umi-length ${umi_length} \\
        --synthetic-length ${synthetic_length} \\
        ${coords_arg}
    
    # Create empty coords file if not provided (to satisfy optional output)
    if [ ! -f "${prefix}_coords_synthetic.csv" ]; then
        touch ${prefix}_coords_synthetic.csv
    fi
    
    cat <<-END_VERSIONS > versions.yml
\t"${task.process}":
\t    python: \$(python3 --version | cut -d' ' -f2)
\tEND_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_synthetic.fastq.gz
    touch ${prefix}_R2_remapped.fastq.gz
    touch ${prefix}_whitelist_synthetic.txt
    touch ${prefix}_barcode_mapping.tsv
    touch ${prefix}_coords_synthetic.csv
    
    cat <<-END_VERSIONS > versions.yml
\t"${task.process}":
\t    python: \$(python3 --version | cut -d' ' -f2)
\tEND_VERSIONS
    """
}
