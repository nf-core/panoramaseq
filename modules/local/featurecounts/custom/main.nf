//custom_featurecounts

process FEATURECOUNTS_CUSTOM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2' :
        'biocontainers/subread:2.0.6--he4a0461_2' }"

    //
    // Now “bam” is a single file, not a list:
    //
    input:
    tuple val(meta), path(bam), path(annotation)

    //
    // “-R BAM” outputs "<prefix>.featureCounts.bam", so we emit that file.
    //
    output:
        tuple val(meta), path("*.featureCounts.tsv"), emit: counts
        tuple val(meta), path("*.featureCounts.tsv.summary"), emit: summary
        tuple val(meta), path("*.featureCounts.bam"), emit: annotated_bam
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def feature_type = task.ext.feature_type ?: 'exon'
    def attribute_type = task.ext.attribute_type ?: 'gene_id'
    def strandedness = task.ext.strandedness ?: '0'  // 0=unstranded, 1=stranded, 2=reverse-stranded
    def paired_end = meta.single_end ? '' : '-p'  // Add -p flag for paired-end data
    """
    featureCounts \\
        -s ${strandedness} \\
        ${paired_end} \\
        -T ${task.cpus} \\
        -t ${feature_type} \\
        -g ${attribute_type} \\
        -a ${annotation} \\
        -o ${prefix}.featureCounts.tsv \\
        -R BAM \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    touch ${prefix}.featureCounts.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
END_VERSIONS
    """
}
