process REORDER_R1_FOR_STARSOLO {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
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
print("Using --soloCBtype String for full barcode support")
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
        
        # Validate read length: need at least 10bp to extract the UMI.
        # Barcode deletions can shorten R1 below 46bp; the barcode is in the
        # QUIK header so only the UMI (first 10bp) is required from the sequence.
        if len(seq) < 10:
            print(f"Warning: Read {records_processed} length {len(seq)} < 10bp, skipping", file=sys.stderr)
            continue
        
        # Extract UMI (first 10bp) and Barcode (remaining chars, variable length
        # when barcode deletions are present).
        umi = seq[:10]
        barcode = seq[10:]   # may be shorter than 36bp due to barcode indels
        
        # Reorder: Barcode first (variable), then UMI (10bp)
        new_seq = barcode + umi
        
        # Quality scores follow same order
        umi_qual = qual[:10]
        bc_qual = qual[10:]
        new_qual = bc_qual + umi_qual
        
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
