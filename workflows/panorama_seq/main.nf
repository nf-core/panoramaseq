/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~    // 11. Count features from BAM using annotation with FEATURECOUNTS_CUSTOM
    //     Uses the provided GTF annotation file
    //     Update metadata to indicate single-end since we only aligned R2
    //     Prepares input as tuple of meta, bam, and annotation file
    custom_featurecounts_input = STAR_ALIGN_LOCAL.out.bam.combine(gtf_file).map { meta, bam, gtf ->
        def new_meta = meta.clone()
        new_meta.single_end = true  // Update to single_end since we only aligned R2
        [new_meta, bam, gtf]
    }~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT2 } from '../../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_panoramaseq_pipeline'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'
include { SAMTOOLS_INDEX as index1; SAMTOOLS_INDEX as index2 } from '../../modules/nf-core/samtools/index/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUTADAPT_PANORAMA as CUTADAPT } from '../../modules/local/cutadapt_panorama/main'
include { Decode_batch } from '../../modules/local/Decoding/main' 
include { CUTADAPT_ADV_PIPE } from '../../modules/local/cutadapt_adv_pipe/main'
include { STAR_ALIGN_LOCAL } from '../../modules/local/star_align/main'
include { FEATURECOUNTS_CUSTOM } from '../../modules/local/FeatCounts/main'   
include { UMI_count } from '../../modules/local/UMI_count/main'
include { SAMTOOLS_SORT_LOCAL } from '../../modules/local/SAMTOOLS_SORT_LOCAL/main'
include { tsv_to_h5ad } from '../../modules/local/anndata/make_h5ad/main' 
include { tsv_to_h5ad_single } from '../../modules/local/anndata/make_h5ad_single/main'
include { CHECK_H5AD } from '../../modules/local/anndata/check_h5ad/main' 

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

workflow panorama_seq {
    
    take:
    valid_data    // channel: samplesheet read in from checked samplesheet process
    star_index    // path: STAR genome index directory
    gtf_file      // path: GTF annotation file
    
    main:

    ch_versions = Channel.empty()
    
    // 1. FastQC on raw input reads (valid_data stage)
    fastqc_raw_input = valid_data.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.id = "raw_${meta.id}"
        [new_meta, reads]
    }
    FASTQC(fastqc_raw_input)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    // 2. Subsample reads for samples with sample_size set in meta using SEQTK_SAMPLE
    seqtk_input = valid_data.filter{ meta, reads -> meta.sample_size != null }
    SEQTK_SAMPLE(seqtk_input)
    ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

    // 3. For samples without sample_size, pass through directly
    passthrough = valid_data.filter{ meta, reads -> meta.sample_size == null }

    // 4. Combine both channels for downstream processing
    // Use concat instead of mix to avoid multi-channel operator issues
    umitools_input = SEQTK_SAMPLE.out.reads.concat(passthrough)

    // 5. Extract UMIs using UMITOOLS_EXTRACT
    umi_extract = UMITOOLS_EXTRACT(umitools_input)
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)

    // 6. Trim reads after UMI extraction using CUTADAPT
    cutadapt_results = CUTADAPT(UMITOOLS_EXTRACT.out.reads)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // 6a. FastQC on reads after first trimming (cutadapt_results stage)
    fastqc_cutadapt_input = CUTADAPT.out.reads.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.id = "cutadapt_${meta.id}"
        [new_meta, reads]
    }
    FASTQC_CUTADAPT(fastqc_cutadapt_input)
    ch_versions = ch_versions.mix(FASTQC_CUTADAPT.out.versions)

    // 7. Decode barcodes using Decode_batch
    // Extract barcode file once (all samples should use the same barcode file)
    ch_barcode_file = CUTADAPT.out.reads.map { meta, reads -> file(meta.barcode_file) }.first()
    decode_results = Decode_batch(
        CUTADAPT.out.reads,
        ch_barcode_file
    )
    ch_versions = ch_versions.mix(Decode_batch.out.versions)

    // 8. Advanced trimming on R2 using CUTADAPT_ADV_PIPE
    cutadapt2_results = CUTADAPT_ADV_PIPE(Decode_batch.out.reads)
    ch_versions = ch_versions.mix(CUTADAPT_ADV_PIPE.out.versions)

    // 8a. FastQC on reads after advanced trimming (cutadapt2_results stage)
    fastqc_cutadapt2_input = CUTADAPT_ADV_PIPE.out.reads.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.id = "cutadapt2_${meta.id}"
        [new_meta, reads]
    }
    FASTQC_CUTADAPT2(fastqc_cutadapt2_input)
    ch_versions = ch_versions.mix(FASTQC_CUTADAPT2.out.versions)

    // 9. Align single-end trimmed R2 fastq using STAR_ALIGN_LOCAL
    //     Uses the provided STAR index directory
    //     Combine the reads with the star_index channel
    star_local_input = CUTADAPT_ADV_PIPE.out.reads.combine(star_index).map { meta, r2fastq, index_dir ->
        [meta, file(r2fastq), index_dir]
    }
    STAR_ALIGN_LOCAL(star_local_input)
    ch_versions = ch_versions.mix(STAR_ALIGN_LOCAL.out.versions)

    // 10. Index the sorted BAM output from STAR_ALIGN_LOCAL using samtools index (index1)
    samtools_index_input = STAR_ALIGN_LOCAL.out.bam
    samtools_index_results = samtools_index_input | index1
    ch_versions = ch_versions.mix(index1.out.versions)

    // 11. Count features from BAM using annotation with FEATURECOUNTS_CUSTOM
    //     Uses the provided GTF annotation file
    //     First update metadata to reflect single-end nature after R2-only alignment
    star_bam_corrected_meta = STAR_ALIGN_LOCAL.out.bam.map { meta, bam ->
        def new_meta = meta.clone()
        new_meta.single_end = true  // Update to single_end since we only aligned R2
        tuple(new_meta, bam)
    }
    //     Then prepare input as tuple of meta, bam, and annotation file
    custom_featurecounts_input = star_bam_corrected_meta.combine(gtf_file).map { meta, bam, gtf ->
        tuple(meta, bam, gtf)
    }
    FEATURECOUNTS_CUSTOM(custom_featurecounts_input)
    ch_versions = ch_versions.mix(FEATURECOUNTS_CUSTOM.out.versions)

    // 12. Sort BAM files after feature counting using SAMTOOLS_SORT_LOCAL
    SAMTOOLS_SORT_LOCAL(FEATURECOUNTS_CUSTOM.out.annotated_bam.map { meta, bam -> tuple(meta, bam) })
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_LOCAL.out.versions)

    // 13. Index the sorted BAM output from SAMTOOLS_SORT_LOCAL using samtools index (index2)
    samtools_index2_input = SAMTOOLS_SORT_LOCAL.out.bam.map { meta, bam -> tuple(meta, bam) }
    samtools_index2_results = samtools_index2_input | index2
    ch_versions = ch_versions.mix(index2.out.versions)

    // 14. Count UMIs using UMI_count
    //     Joins sorted BAM and its index, then passes as tuple to UMI_count
    UMI_count_input = SAMTOOLS_SORT_LOCAL.out.bam.join(index2.out.bai)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }
    UMI_count(UMI_count_input)
    ch_versions = ch_versions.mix(UMI_count.out.versions)

    // 15. Create H5AD files from count TSV files
    if (params.mergecounts) {
        // Merge all count TSV files into single H5AD with spatial coordinates
        ch_counts_with_meta = UMI_count.out.umi_counts
        
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
        
        tsv_to_h5ad(ch_merge_input)
        ch_versions = ch_versions.mix(tsv_to_h5ad.out.versions)
        
        // Collect the merged H5AD for validation
        ch_h5ad_files = tsv_to_h5ad.out.h5ad
        
    } else {
        // Create individual H5AD files for each sample
        ch_single_input = UMI_count.out.umi_counts
            .map { meta, tsv -> 
                def coords_file = file(meta.barcode_file)
                [meta, tsv, coords_file]
            }
        
        tsv_to_h5ad_single(ch_single_input)
        ch_versions = ch_versions.mix(tsv_to_h5ad_single.out.versions)
        
        // Collect individual H5AD files for validation
        ch_h5ad_files = tsv_to_h5ad_single.out.h5ad
    }

    // Optional: Validate H5AD file structure
    if (params.validate_h5ad) {
        CHECK_H5AD(ch_h5ad_files)
        ch_versions = ch_versions.mix(CHECK_H5AD.out.versions)
    }

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
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CUTADAPT.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CUTADAPT2.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.log.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_LOCAL.out.log_final.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_LOCAL.out.log_out.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_CUSTOM.out.summary.collect{it[1]}.ifEmpty([]))

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
