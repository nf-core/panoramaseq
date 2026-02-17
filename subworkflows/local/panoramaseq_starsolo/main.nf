/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUIK_STARSOLO } from '../../../modules/local/quik/starsolo/main'
include { REORDER_R1_FOR_STARSOLO } from '../../../modules/local/reorder_r1/main'
include { REMAP_BARCODES_FOR_STARSOLO } from '../../../modules/local/remap_barcodes/main'
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
    
    // 4. Remap 36bp barcodes to synthetic ≤31bp barcodes for STARsolo compatibility
    //    This preserves full barcode information while bypassing STARsolo's 31bp limit
    REMAP_BARCODES_FOR_STARSOLO(
        REORDER_R1_FOR_STARSOLO.out.reads,
        QUIK_STARSOLO.out.whitelist
    )
    ch_versions = ch_versions.mix(REMAP_BARCODES_FOR_STARSOLO.out.versions.first())
    
    // 5. Combine R2 (cDNA) with remapped R1 (synthetic BC + UMI)
    //    IMPORTANT: STARSOLO module expects [R1, R2] order (it reverses them for STAR)
    ch_starsolo_input = QUIK_STARSOLO.out.r2
        .join(REMAP_BARCODES_FOR_STARSOLO.out.reads)
        .map { meta, r2, r1_synthetic ->
            def new_meta = meta + [solotype: 'CB_UMI_Simple']
            tuple(new_meta, 'CB_UMI_Simple', [r1_synthetic, r2])  // meta, solotype, reads [R1_synthetic, R2]
        }
    
    // 6. Prepare synthetic whitelist channel (STARsolo needs it as a file input)
    ch_whitelist_synthetic = REMAP_BARCODES_FOR_STARSOLO.out.whitelist
    
    // 7. Prepare STAR index channel
    ch_star_index_tuple = ch_star_index
        .map { index -> [[id: 'star_index'], index] }
    
    // 8. Run STARsolo (align + demux + quantify) using synthetic barcodes
    STARSOLO(
        ch_starsolo_input,
        ch_whitelist_synthetic,
        ch_star_index_tuple
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())
    
    // 9. Convert STARsolo output to H5AD with spatial coordinates
    STARSOLO_TO_H5AD(
        STARSOLO.out.counts,      // STARsolo counts directory (*.Solo.out)
        ch_barcode_file
    )
    ch_versions = ch_versions.mix(STARSOLO_TO_H5AD.out.versions.first())
    
    emit:
    h5ad              = STARSOLO_TO_H5AD.out.h5ad                    // channel: [ val(meta), path(h5ad) ]
    star_log_final    = STARSOLO.out.log_final                      // channel: [ val(meta), path(log) ]
    star_log_out      = STARSOLO.out.log_out                        // channel: [ val(meta), path(log) ]
    star_log_progress = STARSOLO.out.log_progress                   // channel: [ val(meta), path(log) ]
    star_counts       = STARSOLO.out.counts                         // channel: [ val(meta), path(Solo.out) ]
    star_summary      = STARSOLO.out.summary                        // channel: [ val(meta), path(Summary.csv) ]
    quik_stats        = QUIK_STARSOLO.out.stats                     // channel: [ val(meta), path(stats) ]
    whitelist         = QUIK_STARSOLO.out.whitelist                 // channel: [ val(meta), path(whitelist) ] - original 36bp
    whitelist_synthetic = REMAP_BARCODES_FOR_STARSOLO.out.whitelist // channel: path(whitelist_synthetic) - synthetic ≤31bp
    barcode_mapping   = REMAP_BARCODES_FOR_STARSOLO.out.mapping     // channel: path(mapping.tsv) - original↔synthetic
    fastqc_zip        = FASTQC.out.zip                              // channel: [ val(meta), path(zip) ]
    versions          = ch_versions                                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
