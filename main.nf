#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================================================
// 0) Import the CHECK_FASTQS process from your fixed module
// ==========================================================================

include { PANORAMASEQ } from './workflows/panoramaseq' 
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome/main'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_panoramaseq_pipeline/main'
include { PANORAMASEQ_COMPLETION } from './subworkflows/local/utils_nfcore_panoramaseq_completion/main'


// ==========================================================================
// 1) Read the sample sheet and build “data” as a flat 3‐element tuple
//    ( meta_map, R1_path, R2_path )
// ==========================================================================


 workflow NFCORE_PANORAMASEQ {  
    
    take:
    valid_data 
 
    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Prepare genome files and generate STAR index
    //
    // Use params.fasta and params.star_gtf if provided to build STAR index
    // Otherwise, rely on params.star_genome_dir (backward compatibility)
    //
    if (params.fasta && params.star_gtf && !params.star_genome_dir) {
        // Generate STAR index from FASTA and GTF
        fasta_file = file(params.fasta, checkIfExists: true)
        gtf_file = file(params.star_gtf, checkIfExists: true)
        
        STAR_GENOMEGENERATE(
            Channel.value([[:], fasta_file]),
            Channel.value([[:], gtf_file])
        )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        
        // Use the generated STAR index
        ch_star_index = STAR_GENOMEGENERATE.out.index.map { meta, index -> index }
        ch_gtf_file = Channel.value(gtf_file)
        
    } else if (params.fasta && params.star_gtf) {
        fasta_file = file(params.fasta, checkIfExists: true)
        gtf_file = file(params.star_gtf, checkIfExists: true)
        ch_additional_fasta = params.additional_fasta ? file(params.additional_fasta, checkIfExists: true) : null

        PREPARE_GENOME(
            fasta_file,
            gtf_file,
            ch_additional_fasta
        )
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
        
        // Use the generated STAR index
        ch_star_index = PREPARE_GENOME.out.index
        ch_gtf_file = PREPARE_GENOME.out.gtf
        
    } else if (params.star_genome_dir && params.star_gtf) {
        // Use pre-built STAR index (backward compatibility)
        ch_star_index = Channel.value(file(params.star_genome_dir, checkIfExists: true))
        ch_gtf_file = Channel.value(file(params.star_gtf, checkIfExists: true))
        
    } else {
        error "ERROR: Either provide --fasta and --star_gtf to build STAR index, or --star_genome_dir and --star_gtf to use existing index"
    }

    // Run main PANORAMASEQ workflow
    // This will handle all the steps defined in the PANORAMASEQ process including quality control, alignment, counting, etc.
    // The PANORAMASEQ process is defined in the workflows/main.nf file
    PANORAMASEQ(
        valid_data,
        ch_star_index,
        ch_gtf_file
    ) // Pass the samples and genome references

    emit:
    multiqc_report = PANORAMASEQ.out.multiqc_report // channel: /path/to/multi
    versions       = ch_versions.mix(PANORAMASEQ.out.versions)
}




workflow {
    main:
    
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )
    
    // Print the output of PIPELINE_INITIALISATION.out.samplesheet
    PIPELINE_INITIALISATION.out.samplesheet.view { "PIPELINE_INITIALISATION.out.samplesheet: $it" }

    // main PANORAMASEQ workflow
    NFCORE_PANORAMASEQ (
         PIPELINE_INITIALISATION.out.samplesheet
    )
    // SUBWORKFLOW: Pipeline completion tasks
    PANORAMASEQ_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PANORAMASEQ.out.multiqc_report
    )


}
// ==========================================================================