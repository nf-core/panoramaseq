//custom_featurecounts

process FEATURECOUNTS_CUSTOM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2'
        : 'biocontainers/subread:2.0.6--he4a0461_2'}"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    featureCounts \\
        -s 1 \\
        -T ${task.cpus} \\
        -t exon \\
        -g gene_id \\
        -a ${annotation} \\
        -o ${prefix}.featureCounts.tsv \\
        -R BAM \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}
