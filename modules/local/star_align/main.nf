process STAR_ALIGN_LOCAL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../nf-core/star/align/environment.yml"
    container "community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --genomeDir $genome_dir \
        --readFilesIn $fastq2 \
        --outFilterType BySJout \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 200 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --limitOutSJcollapsed 5000000 \
        --limitIObufferSize 200000000 200000000 \
        --outSAMattributes NH HI NM MD \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${prefix}_ \
        --limitBAMsortRAM 2000000000

    STAR --version | head -n 1 | sed 's/STAR_//' > versions.yml
    """
}
