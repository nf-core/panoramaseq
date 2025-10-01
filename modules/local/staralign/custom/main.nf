process STAR_ALIGN_LOCAL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"
    input:
    tuple val(meta), path(fastq2), path(genome_dir)

    output:
        tuple val(meta), path("*.bam"), emit: bam
        tuple val(meta), path("*.bai"), optional: true, emit: bai
        tuple val(meta), path("*.cram"), optional: true, emit: cram
        tuple val(meta), path("*.crai"), optional: true, emit: crai
        tuple val(meta), path("*.csi"), optional: true, emit: csi
        path "*_Log.final.out", emit: log_final
        path "*_Log.out", emit: log_out
        path "*_Log.progress.out", emit: log_progress
        path  "versions.yml", emit: versions

    publishDir "${params.outdir}/star", mode: 'copy', overwrite: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --genomeDir $genome_dir \\
        --readFilesIn $fastq2 \\
        --outFileNamePrefix ${prefix}_ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | head -n 1 | sed 's/STAR_//')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_Aligned.sortedByCoord.out.bam
    touch ${prefix}_Log.final.out
    touch ${prefix}_Log.out
    touch ${prefix}_Log.progress.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | head -n 1 | sed 's/STAR_//')
END_VERSIONS
    """
}
