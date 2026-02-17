/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT2 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_panoramaseq_pipeline'
include { UMITOOLS_EXTRACT } from '../modules/nf-core/umitools/extract'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
// include { SAMTOOLS_SORT as SAMTOOLS_SORT_LOCAL} from '../modules/nf-core/samtools/sort/main' //just for testing the nf-core module
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'
include { SAMTOOLS_INDEX as index1; SAMTOOLS_INDEX as index2 } from '../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUTADAPT } from '../modules/nf-core/cutadapt/main'
include { QUIK_BARCODE_CALLING } from '../modules/local/quik/main'
include { CUTADAPT_ADV_PIPE } from '../modules/local/cutadapt_adv_pipe/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { FEATURECOUNTS_CUSTOM } from '../modules/local/featurecounts/custom/main'
include { UMICOUNT } from '../modules/local/umicount/custom/main'
// include { SAMTOOLS_SORT_LOCAL } from '../modules/local/samtoolssort/custom/main'
include { ANNDATA_MAKEH5AD } from '../modules/local/anndata/makeh5ad/main'
include { ANNDATA_MAKEH5AD_SINGLE } from '../modules/local/anndata/makeh5adsingle/main'
include { ANNDATA_CHECKH5AD } from '../modules/local/anndata/checkh5ad/main'
include { PANORAMASEQ_STARSOLO } from '../subworkflows/local/panoramaseq_starsolo/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/nf-core-panoramaseq_logo_light.png", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANORAMASEQ {

    take:
    valid_data    // channel: samplesheet read in from checked samplesheet process
    star_index    // path: STAR genome index directory
    gtf_file      // path: GTF annotation file

    main:

    ch_versions = Channel.empty()

    // 1. FastQC on raw input reads (valid_data stage)
    fastqc_raw_input = valid_data.map { meta, reads ->
        def new_meta = meta + [id:"raw_${meta.id}"]
        [new_meta, reads]
    }
    FASTQC(fastqc_raw_input)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // 2. Subsample reads using SEQTK_SAMPLE if params.sample_size is set
    if (params.sample_size) {
        seqtk_input = valid_data.map { meta, reads ->
            [meta, reads, params.sample_size]
        }
        SEQTK_SAMPLE(seqtk_input)
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())
        ch_workflow_input = SEQTK_SAMPLE.out.reads
    } else {
        // Skip subsampling if sample_size not specified
        ch_workflow_input = valid_data
    }

    //
    // WORKFLOW BRANCHING: Run either BASIC or STARSOLO workflow
    //
    ch_barcode_file = ch_workflow_input.map { meta, reads -> file(meta.barcode_file) }.first()
    
    if (params.workflow_type == 'starsolo') {
        //
        // STARsolo Workflow: QUIK → Reorder R1 → STARsolo → H5AD
        //
        PANORAMASEQ_STARSOLO(
            ch_workflow_input,
            star_index,
            gtf_file,
            ch_barcode_file
        )
        ch_versions = ch_versions.mix(PANORAMASEQ_STARSOLO.out.versions)
        
        // Set outputs for STARsolo workflow
        ch_h5ad_files = PANORAMASEQ_STARSOLO.out.h5ad
        ch_star_logs_final = PANORAMASEQ_STARSOLO.out.star_log_final
        ch_star_logs_out = PANORAMASEQ_STARSOLO.out.star_log_out
        ch_quik_stats = PANORAMASEQ_STARSOLO.out.quik_stats
        ch_fastqc_raw = PANORAMASEQ_STARSOLO.out.fastqc_zip
        
        // Empty channels for unused BASIC workflow outputs
        ch_umi_logs = Channel.empty()
        ch_featurecounts_summary = Channel.empty()
        ch_fastqc_cutadapt = Channel.empty()
        ch_fastqc_cutadapt2 = Channel.empty()
        
    } else {
        //
        // BASIC Workflow: QUIK → UMItools → Cutadapt → STAR → FeatureCounts → UMItools Count → H5AD
        //

    // 3. Decode barcodes using QUIK_BARCODE_CALLING (GPU-accelerated)
    decode_results = QUIK_BARCODE_CALLING(
        ch_workflow_input,
        ch_barcode_file
    )
    ch_versions = ch_versions.mix(QUIK_BARCODE_CALLING.out.versions.first())

    // 4. Extract UMIs using UMITOOLS_EXTRACT
    umi_extract = UMITOOLS_EXTRACT(QUIK_BARCODE_CALLING.out.reads)
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

    // 5. Trim reads after UMI extraction using CUTADAPT
    cutadapt_results = CUTADAPT(UMITOOLS_EXTRACT.out.reads)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    // 6. FastQC on reads after first trimming (cutadapt_results stage)
    fastqc_cutadapt_input = CUTADAPT.out.reads.map { meta, reads ->
        def new_meta = meta + [id: "cutadapt_${meta.id}"]
        [new_meta, reads]
    }
    FASTQC_CUTADAPT(fastqc_cutadapt_input)
    ch_versions = ch_versions.mix(FASTQC_CUTADAPT.out.versions.first())

    // 7. Advanced trimming on R2 using CUTADAPT_ADV_PIPE
    cutadapt2_results = CUTADAPT_ADV_PIPE(CUTADAPT.out.reads)
    ch_versions = ch_versions.mix(CUTADAPT_ADV_PIPE.out.versions.first())

    // 8. FastQC on reads after advanced trimming (cutadapt2_results stage)
    fastqc_cutadapt2_input = CUTADAPT_ADV_PIPE.out.reads.map { meta, reads ->
        def new_meta = meta + [id:"cutadapt2_${meta.id}"]
        [new_meta, reads]
    }
    FASTQC_CUTADAPT2(fastqc_cutadapt2_input)
    ch_versions = ch_versions.mix(FASTQC_CUTADAPT2.out.versions.first())

    // 9. Align single-end trimmed R2 fastq using STAR_ALIGN (nf-core)
    //     Uses the provided STAR index directory and GTF file
    //     Prepare input with single_end flag set to true
    star_input = CUTADAPT_ADV_PIPE.out.reads.map { meta, r2fastq ->
        def new_meta = meta + [single_end: true]
        tuple(new_meta, [r2fastq])  // Wrap reads in list for nf-core module
    }
    STAR_ALIGN(
        star_input,                    // tuple val(meta), path(reads)
        star_index.map { [[:], it] },  // tuple val(meta2), path(index)
        gtf_file.map { [[:], it] },    // tuple val(meta3), path(gtf)
        false,                         // star_ignore_sjdbgtf
        '',                            // seq_platform
        ''                             // seq_center
    )
    // Note: STAR_ALIGN versions are collected automatically via topic emissions

    // 10. Index the sorted BAM output from STAR_ALIGN using samtools index (index1)
    //     Use bam_sorted_aligned output which contains BAM SortedByCoordinate
    samtools_index_input = STAR_ALIGN.out.bam_sorted_aligned
    index1(samtools_index_input)
    ch_versions = ch_versions.mix(index1.out.versions.first())

    // 11. Join BAM and BAI files for featureCounts
    //     featureCounts needs both BAM and index staged together
    //     The join operation matches channels by meta, staging both files in the work directory
    bam_with_index = STAR_ALIGN.out.bam_sorted_aligned
        .join(index1.out.bai, by: 0)  // Join by meta (first element)
        .map { meta, bam, bai ->
            def new_meta = meta + [single_end: true]  // Update to single_end since we only aligned R2
            tuple(new_meta, bam, bai)  // Pass both BAM and BAI
        }

    // 12. Prepare input for FEATURECOUNTS_CUSTOM
    //     featureCounts will have both BAM and BAI staged in its work directory
    //     The BAI is automatically found by featureCounts when it looks for bam_file.bai
    custom_featurecounts_input = bam_with_index
        .combine(gtf_file)
        .map { meta, bam, bai, gtf ->
            tuple(meta, bam, gtf)  // featureCounts input expects (meta, bam, gtf)
        }
    FEATURECOUNTS_CUSTOM(custom_featurecounts_input)
    ch_versions = ch_versions.mix(FEATURECOUNTS_CUSTOM.out.versions.first())

    // 13. Sort BAM files after feature counting using SAMTOOLS_SORT (nf-core module)
    SAMTOOLS_SORT(
        FEATURECOUNTS_CUSTOM.out.annotated_bam.map { meta, bam -> tuple(meta, bam) },
        [[],[]]  // Empty tuple for reference FASTA (not needed)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // 14. Index the sorted BAM output from SAMTOOLS_SORT using samtools index (index2)
    samtools_index2_input = SAMTOOLS_SORT.out.bam.map { meta, bam -> tuple(meta, bam) }
    index2(samtools_index2_input)
    ch_versions = ch_versions.mix(index2.out.versions.first())

    // 15. Count UMIs using UMICOUNT
    //     Joins sorted BAM and its index, then passes as tuple to UMICOUNT
    UMICOUNT_input = SAMTOOLS_SORT.out.bam.join(index2.out.bai)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }
    UMICOUNT(UMICOUNT_input)
    ch_versions = ch_versions.mix(UMICOUNT.out.versions.first())

    // 15. Create H5AD files from count TSV files
    if (params.mergecounts) {
        // Merge all count TSV files into single H5AD with spatial coordinates
        ch_counts_with_meta = UMICOUNT.out.umi_counts

        // Create the merged input by collecting all data and creating a single emission
        ch_merge_input = ch_counts_with_meta
            .collect { meta, tsv -> [meta, tsv] }
            .map { items ->
                // items is a flat list: [meta1, tsv1, meta2, tsv2, ...]
                // So we need to group them back into pairs
                def grouped_items = []
                for (int i = 0; i < items.size(); i += 2) {
                    grouped_items.add([items[i], items[i+1]])  // [meta, tsv]
                }

                def merged_meta = [id: 'merged_counts']
                def tsvs = grouped_items.collect { it[1] }  // Extract all TSV files
                def coords_file = file(items[0].barcode_file)  // Use barcode_file from first sample
                tuple(merged_meta, tsvs, coords_file)
            }

        ANNDATA_MAKEH5AD(ch_merge_input)
        ch_versions = ch_versions.mix(ANNDATA_MAKEH5AD.out.versions.first())

        // Collect the merged H5AD for validation
        ch_h5ad_files = ANNDATA_MAKEH5AD.out.h5ad

    } else {
        // Create individual H5AD files for each sample
        ch_single_input = UMICOUNT.out.umi_counts
            .map { meta, tsv ->
                def coords_file = file(meta.barcode_file)
                [meta, tsv, coords_file]
            }

        ANNDATA_MAKEH5AD_single(ch_single_input)
        ch_versions = ch_versions.mix(ANNDATA_MAKEH5AD_single.out.versions.first())

        // Collect individual H5AD files for validation
        ch_h5ad_files = ANNDATA_MAKEH5AD_single.out.h5ad
    }

    // Optional: Validate H5AD file structure
    if (params.validate_h5ad) {
        ANNDATA_CHECKH5AD(ch_h5ad_files)
        ch_versions = ch_versions.mix(ANNDATA_CHECKH5AD.out.versions.first())
    }

        // Set outputs for BASIC workflow
        ch_star_logs_final = STAR_ALIGN_LOCAL.out.log_final
        ch_star_logs_out = STAR_ALIGN_LOCAL.out.log_out
        ch_umi_logs = UMITOOLS_EXTRACT.out.log
        ch_featurecounts_summary = FEATURECOUNTS_CUSTOM.out.summary
        ch_fastqc_raw = FASTQC.out.zip
        ch_fastqc_cutadapt = FASTQC_CUTADAPT.out.zip
        ch_fastqc_cutadapt2 = FASTQC_CUTADAPT2.out.zip
        ch_quik_stats = QUIK_BARCODE_CALLING.out.stats
    
    }  // End of BASIC vs STARSOLO workflow branching

    // Collect and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_panoramaseq_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // MODULE: MultiQC
    if (!params.skip_multiqc) {
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description  = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_raw.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_cutadapt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_cutadapt2.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_umi_logs.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_logs_final.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_logs_out.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_featurecounts_summary.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        multiqc_report = MULTIQC.out.report.toList()
    } else {
        multiqc_report = Channel.empty()
    }

    emit:
        multiqc_report = multiqc_report // channel: /path/to/multiqc_report.html
        versions       = ch_versions    // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
