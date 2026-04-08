process COLUMBA_RESCUE_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    // All three inputs joined on meta before calling (see subworkflow)
    tuple val(meta), path(sam), path(r1_rejected), path(r2_rejected)

    output:
    tuple val(meta), path("*_R1_rescued.fastq.gz"),    emit: r1_rescued
    tuple val(meta), path("*_R2_rescued.fastq.gz"),    emit: r2_rescued
    tuple val(meta), path("*_rescued_whitelist.txt"),  emit: whitelist
    tuple val(meta), path("*_columba_rescue_stats.txt"), emit: stats
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - <<'PYEOF'
import gzip
import sys

prefix = "${prefix}"
sam_file    = "${sam}"
r1_rejected = "${r1_rejected}"
r2_rejected = "${r2_rejected}"

# ------------------------------------------------------------------
# 1. Parse Columba SAM: collect (read_id -> barcode) for mapped reads.
#    RNAME equals the barcode sequence because we built the index with
#    barcode sequences as FASTA headers.
# ------------------------------------------------------------------
rescued = {}   # read_id -> barcode_sequence
total_in_sam = 0
with open(sam_file) as f:
    for line in f:
        if line.startswith('@'):
            continue
        fields = line.split('\\t')
        if len(fields) < 6:
            continue
        total_in_sam += 1
        read_id = fields[0]
        flag    = int(fields[1])
        rname   = fields[2]
        # Skip unmapped reads (FLAG bit 0x4)
        if flag & 4:
            continue
        # Avoid duplicates: keep first (best) alignment per read
        if read_id not in rescued:
            rescued[read_id] = rname

n_rescued_reads   = len(rescued)
n_unique_barcodes = len(set(rescued.values()))
print(f"SAM records: {total_in_sam}, "
      f"rescued reads: {n_rescued_reads}, "
      f"unique rescued barcodes: {n_unique_barcodes}", flush=True)

# ------------------------------------------------------------------
# 2. Extract rescued R1/R2 reads from the rejected FASTQ and annotate
#    headers with corrected barcodes in QUIK format for downstream
#    barcode remapping compatibility.
# ------------------------------------------------------------------
def extract_reads_with_barcodes(in_path, out_path, rescued_dict):
    """
    Extract reads that were rescued by Columba and update their headers
    to include the corrected barcode in QUIK format:
    @read_id_calledidx_0_<CORRECTED_BARCODE>
    
    This ensures the barcode remapping script can extract the corrected
    barcode from the header and successfully map it.
    """
    n = 0
    opener = gzip.open if str(in_path).endswith('.gz') else open
    with opener(str(in_path), 'rt') as fin, gzip.open(out_path, 'wt') as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            p = fin.readline()
            q = fin.readline()
            
            # Extract read ID (everything after @ before first space/newline)
            read_id = h[1:].split()[0].rstrip()
            
            if read_id in rescued_dict:
                corrected_barcode = rescued_dict[read_id]
                # Rewrite header in QUIK format with corrected barcode
                # Use calledidx_0 as a marker for Columba-rescued reads
                new_header = f"@{read_id}_calledidx_0_{corrected_barcode}\\n"
                fout.write(new_header + s + p + q)
                n += 1
    return n

n_r1 = extract_reads_with_barcodes(r1_rejected, f"{prefix}_R1_rescued.fastq.gz", rescued)
n_r2 = extract_reads_with_barcodes(r2_rejected, f"{prefix}_R2_rescued.fastq.gz", rescued)
print(f"Written: {n_r1} R1 reads, {n_r2} R2 reads", flush=True)

# ------------------------------------------------------------------
# 3. Write rescued whitelist (unique barcodes that were successfully
#    matched by Columba)
# ------------------------------------------------------------------
with open(f"{prefix}_rescued_whitelist.txt", 'w') as f:
    for bc in sorted(set(rescued.values())):
        f.write(bc + '\\n')

# ------------------------------------------------------------------
# 4. Write rescue statistics
# ------------------------------------------------------------------
with open(f"{prefix}_columba_rescue_stats.txt", 'w') as f:
    f.write("Columba Barcode Rescue Summary\\n")
    f.write("=" * 42 + "\\n")
    f.write(f"Sample: {prefix}\\n")
    f.write(f"Total rejected reads submitted to Columba: {total_in_sam}\\n")
    f.write(f"Reads successfully rescued: {n_rescued_reads}\\n")
    f.write(f"Unique rescued barcodes: {n_unique_barcodes}\\n")
    if total_in_sam > 0:
        pct = 100.0 * n_rescued_reads / total_in_sam
        f.write(f"Rescue rate: {pct:.2f}%\\n")

with open("versions.yml", 'w') as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1_rescued.fastq.gz
    touch ${prefix}_R2_rescued.fastq.gz
    touch ${prefix}_rescued_whitelist.txt
    touch ${prefix}_columba_rescue_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //; s/ .*//')
    END_VERSIONS
    """
}
