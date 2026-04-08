process BARCODE_TO_FASTA {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"
    
    input:
    tuple val(meta), path(barcode_coords_csv)
    
    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'PYEOF'
import sys

# Convert barcodes_coords.csv to FASTA format for Columba.
# Use the barcode sequence itself as the FASTA header so that
# the Columba SAM RNAME field directly contains the barcode sequence.
seen = set()
n = 0
with open("${barcode_coords_csv}", 'r') as infile, \\
     open("${prefix}_barcodes.fasta", 'w') as outfile:
    next(infile)  # skip header
    for line in infile:
        barcode = line.strip().split(',')[0]
        if barcode and barcode not in seen:
            outfile.write(f">{barcode}\\n{barcode}\\n")
            seen.add(barcode)
            n += 1

print(f"Wrote {n} unique barcodes to FASTA", flush=True)

with open("versions.yml", 'w') as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
PYEOF
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_barcodes.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
