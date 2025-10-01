//decode Using 2 gpus and high resolution array
process DECODE_BATCH {
    tag "$meta.id" 
    label "gpu_process"
    conda "${moduleDir}/environment.yml"
    container 'docker://pytorch/pytorch:2.7.1-cuda11.8-cudnn9-runtime'
    
    input:
    tuple val(meta), path(reads)
    path barcode_file
    
    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    // Defensive debug print and assertion
    println "[DECODE_BATCH] meta.id: ${meta.id}, N_barcodes: ${meta.N_barcodes}, sample_size: ${meta.sample_size}"
    assert meta.N_barcodes != null : "N_barcodes is null for sample ${meta.id}"
    assert meta.N_barcodes > 0 : "N_barcodes is not positive for sample ${meta.id}"

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gpus = params.gpus ?: 2
    def cmds = []
    for (int i = 1; i <= gpus; i++) {
        int cuda = (i-1) % gpus
        cmds << "Demux2_batch.py -R1 ${reads[0]} -R2 ${reads[1]} -N ${meta.N_barcodes} -M ${meta.len_barcode} -B ${barcode_file} --Ntriage ${meta.Ntriage} --Nthresh ${meta.Nthresh} --out-prefix ${meta.id} -C ${cuda} --Seg_id ${i}/${gpus} &"
        cmds << "pid${i-1}=\$!"
    }
    def waits = (0..<gpus).collect { "wait \$pid${it}" }.join("\n")
    def catR1 = "cat *${gpus}_R1.fastq.gz > ${prefix}.R1.fastq.gz"
    def catR2 = "cat *${gpus}_R2.fastq.gz > ${prefix}.R2.fastq.gz"
    def rm = "rm -f *${gpus}_R1.fastq.gz *${gpus}_R2.fastq.gz"
    """
    echo CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}
    python -c "import torch; print('torch.cuda.is_available():', torch.cuda.is_available())"
    ${cmds.join("\n")}
    ${waits}
    ${catR1}
    ${catR2}
    ${rm}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        torch: \$(python -c "import torch; print(torch.__version__)")
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create stub output files for testing
    touch ${prefix}.R1.fastq.gz
    touch ${prefix}.R2.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub-version"
        torch: "stub-version"
END_VERSIONS
    """
}