process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    label "process_single"

    input:
    path samplesheet

    output:
    path "samplesheet_checked.log", emit: log
    path "samplesheet_valid.csv", emit: valid

    script:
    """
    python3 -c "import sys, csv; required = ['sample', 'fastq_1', 'fastq_2', 'N_barcodes', 'barcode_file', 'sample_size', 'len_barcode', 'Nthresh', 'Ntriage'];\nwith open('${samplesheet}') as f:\n    reader = csv.DictReader(f)\n    fieldnames = [fn.strip().replace('\\r','').replace('\\n','') for fn in reader.fieldnames]\n    with open('samplesheet_checked.log', 'w') as log:\n        log.write('Detected columns: ' + str(fieldnames) + '\\n')\n        missing = [col for col in required if col not in fieldnames]\n        if missing:\n            log.write('ERROR: Missing columns: ' + ','.join(missing) + '\\n')\n            sys.exit(1)\n        else:\n            log.write('Samplesheet OK\\n')\n            with open('samplesheet_valid.csv', 'w', newline='') as out:\n                writer = csv.DictWriter(out, fieldnames=required)\n                writer.writeheader()\n                for i, row in enumerate(reader, 1):\n                    if not row.get('Ntriage') or str(row.get('Ntriage')).strip() == '':\n                        log.write(f'ERROR: Row {i} missing Ntriage value\\n')\n                        sys.exit(2)\n                    writer.writerow({col: row.get(col, '') for col in required})"
    """
}
