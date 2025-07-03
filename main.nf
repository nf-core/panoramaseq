#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================================================
// 0) Import the CHECK_FASTQS process from your fixed module
// ==========================================================================

include { panorama_seq } from './workflows/panorama_seq/main' 
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome/main'
include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_panoramaseq_pipeline/main'

log.info """
\033[1;34m _____                                                                                      _____ 
\033[1;34m( ___ )------------------------------------------------------------------------------------( ___ )
\033[1;34m |   |                                                                                      |   | 
\033[1;34m |   |                                                                                      |   | 
\033[1;34m |   |                                                                                      |   | 
\033[1;34m |   |    ______                                                                            |   | 
\033[1;34m |   |   (_____ \\                                                                           |   | 
\033[1;34m |   |    _____) )  ____  ____    ___    ____   ____  ____    ____     ___   ____   ____    |   | 
\033[1;34m |   |   |  ____/  / _  ||  _ \\  / _ \\  / ___) / _  ||    \\  / _  |   /___) / _  ) / _  |   |   | 
\033[1;34m |   |   | |      ( ( | || | | || |_| || |    ( ( | || | | |( ( | |  |___ |( (/ / | | | |   |   | 
\033[1;34m |   |   |_|       \\_||_||_| |_| \\___/ |_|     \\_||_||_|_|_| \\_||_|  (___/  \\____) \\_|| |   |   | 
\033[1;34m |   |                                                                                |_|   |   | 
\033[1;34m |   |                                                                                      |   | 
\033[1;34m |   |                              Spatial Transcriptomics                                 |   | 
\033[1;34m |___|                                                                                      |___| 
\033[1;34m(_____)------------------------------------------------------------------------------------(_____)
                                                                                                                     
   \033[1;34m ================================       
   \033[1;34m             Author:
   \033[1;34m Franco Alexander Poma Soto
   \033[1;34m OncoRNA Lab
   \033[1;34m Ghent University
   \033[1;34m ================================
    """.stripIndent()

// ==========================================================================
// 1) Read the sample sheet and build “data” as a flat 3‐element tuple
//    ( meta_map, R1_path, R2_path )
// ==========================================================================


 workflow NFCORE_PANORAMAseq {  
    
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
    if (params.fasta && params.star_gtf) {
        ch_fasta = file(params.fasta, checkIfExists: true)
        ch_gtf = file(params.star_gtf, checkIfExists: true)
        ch_additional_fasta = params.additional_fasta ? file(params.additional_fasta, checkIfExists: true) : null

        PREPARE_GENOME(
            ch_fasta,
            ch_gtf,
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
    panorama_seq(
        valid_data,
        ch_star_index,
        ch_gtf_file
    ) // Pass the samples and genome references

    emit:
    multiqc_report = panorama_seq.out.multiqc_report // channel: /path/to/multi
    versions       = ch_versions.mix(panorama_seq.out.versions)
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
    NFCORE_PANORAMAseq (
         PIPELINE_INITIALISATION.out.samplesheet
    )
    // SUBWORKFLOW: Pipeline completion tasks
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PANORAMAseq.out.multiqc_report
    )


}
// ==========================================================================