/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUIK_STARSOLO } from '../../../modules/local/quik_starsolo/main'
include { REORDER_R1_FOR_STARSOLO } from '../../../modules/local/reorder_r1/main'
include { REMAP_BARCODES_FOR_STARSOLO } from '../../../modules/local/remap_barcodes/main'
include { STARSOLO } from '../../../modules/nf-core/star/starsolo/main'
include { STARSOLO_TO_H5AD } from '../../../modules/local/starsolo_to_h5ad/main'
include { STARSOLO_HEATMAP } from '../../../modules/local/starsolo_heatmap/main'
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
    
    // 3. Conditionally reorder R1 based on read structure
    //    - UMI_BC (e.g., Visium): swap UMI(10bp) + BC(36bp) → BC(36bp) + UMI(10bp)
    //    - BC_UMI (e.g., custom prep): already in correct order, skip reordering
    if (params.read_structure == 'UMI_BC') {
        REORDER_R1_FOR_STARSOLO(QUIK_STARSOLO.out.r1)
        ch_versions = ch_versions.mix(REORDER_R1_FOR_STARSOLO.out.versions.first())
        ch_r1_ordered = REORDER_R1_FOR_STARSOLO.out.reads
    } else {
        // BC_UMI: already in correct order (BC first, then UMI)
        ch_r1_ordered = QUIK_STARSOLO.out.r1
    }
    
    // 4. Remap 36bp barcodes to synthetic ≤31bp barcodes for STARsolo compatibility
    //    This preserves full barcode information while bypassing STARsolo's 31bp limit
    //    Join ordered R1 with whitelist and R2 by meta.id
    //    ALSO generate synthetic coordinate file to avoid expensive reverse mapping during heatmap
    ch_for_remap = ch_r1_ordered
        .join(QUIK_STARSOLO.out.whitelist)
        .join(QUIK_STARSOLO.out.r2)
    
    REMAP_BARCODES_FOR_STARSOLO(
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, r1] },
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, whitelist] },
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, r2] },
        ch_barcode_file  // Coordinates file to generate synthetic coords
    )
    ch_versions = ch_versions.mix(REMAP_BARCODES_FOR_STARSOLO.out.versions.first())
    
    // 5. Combine filtered R2 (cDNA) with remapped R1 (synthetic BC + UMI)
    //    IMPORTANT: STARSOLO module expects [R1, R2] order (it reverses them for STAR)
    //    Both R1 and R2 are now filtered to have matching read counts
    ch_starsolo_input = REMAP_BARCODES_FOR_STARSOLO.out.reads_r2
        .join(REMAP_BARCODES_FOR_STARSOLO.out.reads)
        .map { meta, r2_remapped, r1_synthetic ->
            def new_meta = meta + [solotype: 'CB_UMI_Simple']
            tuple(new_meta, 'CB_UMI_Simple', [r1_synthetic, r2_remapped])  // meta, solotype, reads [R1_synthetic, R2_remapped]
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
    
    // 10. Generate spatial UMI count heatmap using synthetic barcode coordinates (optimized!)
    //     Using synthetic coords directly avoids expensive reverse mapping
    //     Create a dummy mapping channel (not actually used by the heatmap script)
    ch_heatmap_input = STARSOLO.out.counts
    
    // Use synthetic coords for optimized heatmap generation (no reverse mapping needed)
    ch_synthetic_coords = REMAP_BARCODES_FOR_STARSOLO.out.coords
    
    // Create dummy mapping channel to satisfy process signature (mapping is optional in the script now)
    ch_dummy_mapping = ch_heatmap_input
        .map { meta, counts -> 
            def dummy_file = file("${workflow.workDir}/NO_FILE")
            dummy_file.text = ""
            [meta, dummy_file]
        }
    
    STARSOLO_HEATMAP(
        ch_heatmap_input,
        ch_dummy_mapping,
        ch_synthetic_coords
    )
    ch_versions = ch_versions.mix(STARSOLO_HEATMAP.out.versions.first())
    
    emit:
    h5ad              = STARSOLO_TO_H5AD.out.h5ad                    // channel: [ val(meta), path(h5ad) ]
    heatmap           = STARSOLO_HEATMAP.out.heatmap                 // channel: [ val(meta), path(png) ]
    heatmap_data      = STARSOLO_HEATMAP.out.data                    // channel: [ val(meta), path(tsv) ]
    heatmap_stats     = STARSOLO_HEATMAP.out.stats                   // channel: [ val(meta), path(json) ]
    star_log_final    = STARSOLO.out.log_final                      // channel: [ val(meta), path(log) ]
    star_log_out      = STARSOLO.out.log_out                        // channel: [ val(meta), path(log) ]
    star_log_progress = STARSOLO.out.log_progress                   // channel: [ val(meta), path(log) ]
    star_counts       = STARSOLO.out.counts                         // channel: [ val(meta), path(Solo.out) ]
    star_summary      = STARSOLO.out.summary                        // channel: [ val(meta), path(Summary.csv) ]
    quik_stats        = QUIK_STARSOLO.out.stats                     // channel: [ val(meta), path(stats) ]
    whitelist         = QUIK_STARSOLO.out.whitelist                 // channel: [ val(meta), path(whitelist) ] - original 36bp
    whitelist_synthetic = REMAP_BARCODES_FOR_STARSOLO.out.whitelist // channel: path(whitelist_synthetic) - synthetic ≤31bp
    barcode_mapping   = REMAP_BARCODES_FOR_STARSOLO.out.mapping     // channel: [ val(meta), path(mapping.tsv) ] - original↔synthetic
    coords_synthetic  = REMAP_BARCODES_FOR_STARSOLO.out.coords      // channel: path(coords_synthetic.csv) - synthetic barcode coords
    fastqc_zip        = FASTQC.out.zip                              // channel: [ val(meta), path(zip) ]
    versions          = ch_versions                                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
