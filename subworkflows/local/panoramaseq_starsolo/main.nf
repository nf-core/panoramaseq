/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUIK_STARSOLO } from '../../../modules/local/quik/starsolo/main'
include { REORDER_R1_FOR_STARSOLO } from '../../../modules/local/reorder_r1/main'
include { STARSOLO } from '../../../modules/nf-core/star/starsolo/main'
include { STARSOLO_TO_H5AD } from '../../../modules/local/starsolo_to_h5ad/main'
include { FASTQC } from '../../../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN STARSOLO WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANORAMASEQ_STARSOLO {
    
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_star_index     // channel: /path/to/star/index
    ch_gtf_file       // channel: /path/to/annotation.gtf
    ch_barcode_file   // channel: /path/to/barcodes.csv
    
    main:
    
    ch_versions = Channel.empty()
    
    // 1. FastQC on raw input reads
    fastqc_raw_input = ch_reads.map { meta, reads ->
        def new_meta = meta + [id:"raw_${meta.id}"]
        [new_meta, reads]
    }
    FASTQC(fastqc_raw_input)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    
    // 2. QUIK barcode calling (generates whitelist)
    QUIK_STARSOLO(ch_reads, ch_barcode_file)
    ch_versions = ch_versions.mix(QUIK_STARSOLO.out.versions.first())
    
    // 3. Reorder R1: swap UMI(10bp) + Barcode(36bp) to Barcode(36bp) + UMI(10bp)
    REORDER_R1_FOR_STARSOLO(QUIK_STARSOLO.out.r1)
    ch_versions = ch_versions.mix(REORDER_R1_FOR_STARSOLO.out.versions.first())
    
    // 4. Combine R2 (cDNA) with reordered R1 (BC+UMI)
    //    IMPORTANT: STARsolo expects [R2, R1] order with solotype in meta
    ch_starsolo_input = QUIK_STARSOLO.out.r2
        .join(REORDER_R1_FOR_STARSOLO.out.reads)
        .map { meta, r2, r1 ->
            def new_meta = meta + [solotype: 'CB_UMI_Simple']
            tuple(new_meta, 'CB_UMI_Simple', [r2, r1])  // meta, solotype, reads [R2, R1]
        }
    
    // 5. Prepare whitelist channel (STARsolo needs it as a file input)
    ch_whitelist = QUIK_STARSOLO.out.whitelist
    
    // 6. Prepare STAR index channel
    ch_star_index_tuple = ch_star_index
        .map { index -> [[id: 'star_index'], index] }
    
    // 7. Run STARsolo (align + demux + quantify)
    STARSOLO(
        ch_starsolo_input,
        ch_whitelist,
        ch_star_index_tuple
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())
    
    // 8. Convert STARsolo output to H5AD with spatial coordinates
    STARSOLO_TO_H5AD(
        STARSOLO.out.tab_gene,    // STARsolo Gene output directory
        ch_barcode_file
    )
    ch_versions = ch_versions.mix(STARSOLO_TO_H5AD.out.versions.first())
    
    emit:
    h5ad              = STARSOLO_TO_H5AD.out.h5ad           // channel: [ val(meta), path(h5ad) ]
    star_log_final    = STARSOLO.out.log_final             // channel: [ val(meta), path(log) ]
    star_log_out      = STARSOLO.out.log_out               // channel: [ val(meta), path(log) ]
    star_log_progress = STARSOLO.out.log_progress          // channel: [ val(meta), path(log) ]
    star_bam          = STARSOLO.out.bam                   // channel: [ val(meta), path(bam) ]
    quik_stats        = QUIK_STARSOLO.out.stats            // channel: [ val(meta), path(stats) ]
    whitelist         = QUIK_STARSOLO.out.whitelist        // channel: path(whitelist)
    fastqc_zip        = FASTQC.out.zip                     // channel: [ val(meta), path(zip) ]
    versions          = ch_versions                         // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
