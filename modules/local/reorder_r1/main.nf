process REORDER_R1_FOR_STARSOLO {
    tag "${meta.id}"
    label 'process_medium'
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"
    
    input:
    tuple val(meta), path(r1)
    
    output:
    tuple val(meta), path("*_reordered.fastq.gz"), emit: reads
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 <<'EOF'
import gzip
import sys
import subprocess

print("=== Reordering R1 for STARsolo ===")
print("Input: UMI(10bp) + Barcode(36bp)")
print("Output: Barcode(36bp) + UMI(10bp)")
print("===================================")

records_processed = 0
records_reordered = 0

with gzip.open('${r1}', 'rt') as infile, \\
     gzip.open('${prefix}_R1_reordered.fastq.gz', 'wt') as outfile:
    
    while True:
        # Read 4 lines (one FASTQ record)
        header = infile.readline()
        if not header: 
            break
        
        seq = infile.readline().strip()
        plus = infile.readline()
        qual = infile.readline().strip()
        
        records_processed += 1
        
        # Validate read length (should be at least 46bp: 10bp UMI + 36bp BC)
        if len(seq) < 46:
            print(f"Warning: Read {records_processed} length {len(seq)} < 46bp, skipping", file=sys.stderr)
            continue
        
        # Extract UMI (first 10bp) and Barcode (next 36bp)
        umi = seq[:10]
        barcode = seq[10:46]
        rest = seq[46:]  # Any remaining sequence (usually empty)
        
        # Reorder: Barcode first, then UMI, then rest
        new_seq = barcode + umi + rest
        
        # Quality scores follow same order
        umi_qual = qual[:10]
        bc_qual = qual[10:46]
        rest_qual = qual[46:]
        new_qual = bc_qual + umi_qual + rest_qual
        
        # Write reordered record
        outfile.write(header)
        outfile.write(new_seq + '\\n')
        outfile.write(plus)
        outfile.write(new_qual + '\\n')
        
        records_reordered += 1

print(f"Processed: {records_processed} reads")
print(f"Reordered: {records_reordered} reads")
print(f"Skipped: {records_processed - records_reordered} reads")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_reordered.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
