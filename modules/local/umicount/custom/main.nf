process UMICOUNT {
    tag "$meta.id"
    label "process_long"

    // conda "${moduleDir}/environment.yml"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"
    
    input:
    tuple val(meta), path(bam), path(bam_index)
    
    output:
    tuple val(meta), path("*.tsv.gz"), emit: umi_counts
    path('logs_umi/*.{txt,log}'), emit: log_files
    path "versions.yml", emit: versions

    script:
    """
    mkdir logs_umi
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${bam} -S ${meta.id}_counts.tsv.gz --log=logs_umi/${meta.id}_umiCounts.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$(umi_tools --version | sed '/version:/!d; s/.*: //')
END_VERSIONS
    """

    stub:
    """
    mkdir logs_umi
    # Create a valid gzipped file instead of empty file
    echo -e "gene\\tcell\\tcount\\tumi_count" | gzip > ${meta.id}_counts.tsv.gz
    echo "Stub UMI counting log" > logs_umi/${meta.id}_umiCounts.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: "stub-version"
END_VERSIONS
    """

}