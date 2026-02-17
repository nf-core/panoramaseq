process REMAP_BARCODES_FOR_STARSOLO {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"
    
    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(whitelist)
    
    output:
    tuple val(meta), path("*_synthetic.fastq.gz"), emit: reads
    path "*_whitelist_synthetic.txt"              , emit: whitelist
    path "*_barcode_mapping.tsv"                  , emit: mapping
    path "versions.yml"                           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_start = task.ext.barcode_start ?: 0
    def barcode_length = task.ext.barcode_length ?: 36
    def umi_length = task.ext.umi_length ?: 10
    def synthetic_length = task.ext.synthetic_length ?: 25
    """
    remap_barcodes_for_starsolo.py \\
        --whitelist ${whitelist} \\
        --fastq ${reads} \\
        --output-fastq ${prefix}_synthetic.fastq.gz \\
        --output-whitelist ${prefix}_whitelist_synthetic.txt \\
        --output-mapping ${prefix}_barcode_mapping.tsv \\
        --barcode-start ${barcode_start} \\
        --barcode-length ${barcode_length} \\
        --umi-length ${umi_length} \\
        --synthetic-length ${synthetic_length}
    
    cat <<-END_VERSIONS > versions.yml
\t"${task.process}":
\t    python: \$(python3 --version | cut -d' ' -f2)
\tEND_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_synthetic.fastq.gz
    touch ${prefix}_whitelist_synthetic.txt
    touch ${prefix}_barcode_mapping.tsv
    
    cat <<-END_VERSIONS > versions.yml
\t"${task.process}":
\t    python: \$(python3 --version | cut -d' ' -f2)
\tEND_VERSIONS
    """
}
