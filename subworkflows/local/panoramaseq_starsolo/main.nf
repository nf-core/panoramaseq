/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUIK_STARSOLO }              from '../../../modules/local/quik_starsolo/main'
include { REORDER_R1_FOR_STARSOLO }    from '../../../modules/local/reorder_r1/main'
include { REMAP_BARCODES_FOR_STARSOLO } from '../../../modules/local/remap_barcodes/main'
include { STARSOLO }                   from '../../../modules/nf-core/star/starsolo/main'
include { STARSOLO_TO_H5AD }           from '../../../modules/local/starsolo_to_h5ad/main'
include { STARSOLO_HEATMAP }           from '../../../modules/local/starsolo_heatmap/main'
include { FASTQC }                     from '../../../modules/nf-core/fastqc/main'

// Columba barcode rescue modules (only used when params.enable_columba_rescue = true)
include { COLUMBA_BUILD }              from '../../../modules/local/columba/build/main'
include { BARCODE_TO_FASTA }           from '../../../modules/local/columba/makefasta/main'
include { COLUMBA_INDEX }              from '../../../modules/local/columba/index/main'
include { COLUMBA_ALIGN }              from '../../../modules/local/columba/align/main'
include { COLUMBA_RESCUE_READS }       from '../../../modules/local/columba/rescue_reads/main'
include { MERGE_RESCUED_READS }        from '../../../modules/local/columba/merge_rescued/main'

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

    // 2. QUIK GPU barcode calling (generates whitelist + optionally rejected reads)
    QUIK_STARSOLO(ch_reads, ch_barcode_file)
    ch_versions = ch_versions.mix(QUIK_STARSOLO.out.versions.first())

    // 3. Conditionally reorder R1 based on read structure
    if (params.read_structure == 'UMI_BC') {
        REORDER_R1_FOR_STARSOLO(QUIK_STARSOLO.out.r1)
        ch_versions = ch_versions.mix(REORDER_R1_FOR_STARSOLO.out.versions.first())
        ch_r1_ordered = REORDER_R1_FOR_STARSOLO.out.reads
    } else {
        ch_r1_ordered = QUIK_STARSOLO.out.r1
    }

    // 4. Optional Columba barcode rescue:
    //    Aligns QUIK-rejected reads against the reference barcode set using
    //    approximate string matching, rescuing reads with sequencing errors
    //    that QUIK's threshold could not correct.
    //
    //    Singularity profile: uses pre-built binaries from the container image.
    //    Conda profile: requires params.columba_repo pointing to a local clone of
    //                   the Columba repository (https://github.com/biointec/columba).
    //                   The build script inside the repo will be invoked automatically
    //                   if the binaries have not been compiled yet.
    if (params.enable_columba_rescue) {

        // 4a. Build (or copy from container) Columba binaries — runs once per pipeline.
        //     .first() makes the directory channel reusable across multiple samples.
        ch_columba_repo = Channel.value(params.columba_repo ?: 'NOT_SET')
        COLUMBA_BUILD(ch_columba_repo)
        ch_versions = ch_versions.mix(COLUMBA_BUILD.out.versions.first())
        // Wrap in a value-like channel so it broadcasts to all samples
        ch_binaries_dir = COLUMBA_BUILD.out.binaries_dir.first()

        // 4b. Convert barcode coords CSV → FASTA for indexing (per sample).
        //     The barcode file is shared across all samples (same Visium HD whitelist).
        ch_barcode_for_fasta = ch_reads
            .map    { meta, reads -> meta }
            .combine(ch_barcode_file)
            .map    { meta, barcode_f -> [meta, barcode_f] }
        BARCODE_TO_FASTA(ch_barcode_for_fasta)
        ch_versions = ch_versions.mix(BARCODE_TO_FASTA.out.versions.first())

        // 4c. Build Columba index from barcode FASTA (per sample)
        COLUMBA_INDEX(
            BARCODE_TO_FASTA.out.fasta,
            ch_binaries_dir
        )
        ch_versions = ch_versions.mix(COLUMBA_INDEX.out.versions.first())

        // 4d. Align rejected reads against the barcode index (per sample)
        //     Input channel: [ meta, [ r1_rejected ], index_files ]
        ch_align_input = QUIK_STARSOLO.out.r1_rejected
            .join(COLUMBA_INDEX.out.index, by: 0)
            .map { meta, r1_rej, idx_files -> [meta, [r1_rej], idx_files] }

        COLUMBA_ALIGN(
            ch_align_input,
            ch_binaries_dir
        )
        ch_versions = ch_versions.mix(COLUMBA_ALIGN.out.versions.first())

        // 4e. Parse SAM → rescued R1/R2 reads + rescued whitelist.
        //     Join all three channels on meta to guarantee correct pairing.
        ch_rescue_input = COLUMBA_ALIGN.out.sam
            .join(QUIK_STARSOLO.out.r1_rejected, by: 0)
            .join(QUIK_STARSOLO.out.r2_rejected,  by: 0)
            .map { meta, sam, r1_rej, r2_rej -> [meta, sam, r1_rej, r2_rej] }
        COLUMBA_RESCUE_READS(ch_rescue_input)
        ch_versions = ch_versions.mix(COLUMBA_RESCUE_READS.out.versions.first())

        // 4f. Merge QUIK-filtered + Columba-rescued streams.
        //     Join all six per-sample channels on meta before calling the process.
        ch_merge_input = ch_r1_ordered
            .join(QUIK_STARSOLO.out.r2,                by: 0)
            .join(QUIK_STARSOLO.out.whitelist,          by: 0)
            .join(COLUMBA_RESCUE_READS.out.r1_rescued,  by: 0)
            .join(COLUMBA_RESCUE_READS.out.r2_rescued,  by: 0)
            .join(COLUMBA_RESCUE_READS.out.whitelist,   by: 0)
            .map { meta, r1f, r2f, wl_quik, r1r, r2r, wl_rescued ->
                [meta, r1f, r2f, wl_quik, r1r, r2r, wl_rescued]
            }
        MERGE_RESCUED_READS(ch_merge_input)
        ch_versions = ch_versions.mix(MERGE_RESCUED_READS.out.versions.first())

        ch_r1_final   = MERGE_RESCUED_READS.out.r1
        ch_r2_final   = MERGE_RESCUED_READS.out.r2
        ch_whitelist  = MERGE_RESCUED_READS.out.whitelist

    } else {
        ch_r1_final   = ch_r1_ordered
        ch_r2_final   = QUIK_STARSOLO.out.r2
        ch_whitelist  = QUIK_STARSOLO.out.whitelist
    }

    // 5. Remap 36bp → synthetic ≤31bp barcodes
    ch_for_remap = ch_r1_final
        .join(ch_whitelist)
        .join(ch_r2_final)

    REMAP_BARCODES_FOR_STARSOLO(
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, r1] },
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, whitelist] },
        ch_for_remap.map { meta, r1, whitelist, r2 -> [meta, r2] },
        ch_barcode_file
    )
    ch_versions = ch_versions.mix(REMAP_BARCODES_FOR_STARSOLO.out.versions.first())

    // 6. Prepare STARsolo input with synthetic barcodes
    ch_starsolo_input = REMAP_BARCODES_FOR_STARSOLO.out.reads_r2
        .join(REMAP_BARCODES_FOR_STARSOLO.out.reads)
        .map { meta, r2_remapped, r1_synthetic ->
            def new_meta = meta + [solotype: 'CB_UMI_Simple']
            tuple(new_meta, 'CB_UMI_Simple', [r1_synthetic, r2_remapped])
        }

    ch_whitelist_starsolo = REMAP_BARCODES_FOR_STARSOLO.out.whitelist
    ch_coords_for_heatmap = REMAP_BARCODES_FOR_STARSOLO.out.coords

    // 7. Prepare STAR index channel
    ch_star_index_tuple = ch_star_index
        .map { index -> [[id: 'star_index'], index] }

    // 8. Run STARsolo (align + demux + quantify)
    STARSOLO(
        ch_starsolo_input,
        ch_whitelist_starsolo,
        ch_star_index_tuple
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())

    // 9. Convert STARsolo output to H5AD with spatial coordinates
    STARSOLO_TO_H5AD(
        STARSOLO.out.counts,
        ch_barcode_file  // Always use original coords for H5AD
    )
    ch_versions = ch_versions.mix(STARSOLO_TO_H5AD.out.versions.first())

    // 10. Generate spatial UMI count heatmap
    ch_heatmap_input = STARSOLO.out.counts

    // Create dummy mapping channel (not used by heatmap script)
    ch_dummy_mapping = ch_heatmap_input
        .map { meta, counts ->
            def dummy_file = file("${workflow.workDir}/NO_FILE")
            dummy_file.text = ""
            [meta, dummy_file]
        }

    STARSOLO_HEATMAP(
        ch_heatmap_input,
        ch_dummy_mapping,
        ch_coords_for_heatmap
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
    whitelist         = ch_whitelist                                // channel: [ val(meta), path(whitelist) ]
    fastqc_zip        = FASTQC.out.zip                              // channel: [ val(meta), path(zip) ]
    versions          = ch_versions                                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
