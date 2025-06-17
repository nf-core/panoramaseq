process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
    path  "versions.yml"           , emit: versions

    publishDir "${params.outdir}/star", mode: 'copy', overwrite: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram' ? 'crai' : args.contains('-c') ? 'csi' : 'bai'
    def index_file = file(input).getName() + '.' + extension
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram' ?
                    "crai" : args.contains("-c") ?  "csi" : "bai"
    def index_file = file(input).getName() + '.' + extension
    """
    touch $index_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
