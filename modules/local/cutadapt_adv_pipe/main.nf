// Custom advanced, piped cutadapt for R2 only
process CUTADAPT_ADV_PIPE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/1758869538eb8e658077cc14cd7a4e76fd9b6d73d3a68f85a70bf292e39e27c5/data' :
        'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' }"
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*.trim2.R2.fastq.gz"), emit: reads
    path "versions.yml", emit: versions
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${reads[1]} \
    | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - \
    | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - \
    | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${prefix}.trim2.R2.fastq.gz -
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub output files for testing
    touch ${prefix}.trim2.R2.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: "stub-version"
END_VERSIONS
    """
}
