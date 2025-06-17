#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ==========================================================================
// 0) Import the CHECK_FASTQS process from your fixed module
// ==========================================================================

include { CHECK_SAMPLESHEET } from './modules/local/checksamplesheet'
include { ST_main } from './workflows/ST_main/main' 
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_panoramaseq_pipeline/main'
log.info """
 _____                                                                                      _____ 
( ___ )------------------------------------------------------------------------------------( ___ )
 |   |                                                                                      |   | 
 |   |                                                                                      |   | 
 |   |                                                                                      |   | 
 |   |    ______                                                                            |   | 
 |   |   (_____ \\                                                                           |   | 
 |   |    _____) )  ____  ____    ___    ____   ____  ____    ____     ___   ____   ____    |   | 
 |   |   |  ____/  / _  ||  _ \\  / _ \\  / ___) / _  ||    \\  / _  |   /___) / _  ) / _  |   |   | 
 |   |   | |      ( ( | || | | || |_| || |    ( ( | || | | |( ( | |  |___ |( (/ / | | | |   |   | 
 |   |   |_|       \\_||_||_| |_| \\___/ |_|     \\_||_||_|_|_| \\_||_|  (___/  \\____) \\_|| |   |   | 
 |   |                                                                                |_|   |   | 
 |   |                                                                                      |   | 
 |   |                                                                                      |   | 
 |___|                                                                                      |___| 
(_____)------------------------------------------------------------------------------------(_____)
                                                                                                                    
        LIST OF PARAMETERS
    ================================
                INPUT
    Sample Sheet     : $params.input 
    Results-folder   : $params.outdir
    ================================
    
    ================================
            UMI-tools
    extract method  : $params.umitools_extract_method
    UMI-pattern        : $params.umitools_bc_pattern          
    ================================   
            Decoding
    number of gpus   : $params.gpus
    ================================       
                STAR
    Reference genome : $params.star_genome_dir
    GTF-file         : $params.star_gtf
    ================================
    """.stripIndent()

// ==========================================================================
// 1) Read the sample sheet and build “data” as a flat 3‐element tuple
//    ( meta_map, R1_path, R2_path )
// ==========================================================================


 workflow NFCORE_PANORAMAseq {  
    
    take:
    valid_data 
 
    main:

    // Run main PANORAMASEQ workflow
    // This will handle all the steps defined in the PANORAMASEQ process including quality control, alignment, counting, etc.
    // The PANORAMASEQ process is defined in the workflows/main.nf file
    ST_main(valid_data) // Pass the samples

    // emit:
    // multiqc_report = ST_main.out.multiqc_report // channel: /path/to/multi
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

    // Comment out the main workflow for now
    NFCORE_PANORAMAseq (
         PIPELINE_INITIALISATION.out.samplesheet
    )
}
// ==========================================================================