//
// Prepare genome: process genome files and generate STAR index
//

include { GUNZIP as GUNZIP_FASTA    } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL } from '../../../modules/nf-core/gunzip/main'
include { STAR_GENOMEGENERATE       } from '../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME {

    take:
    fasta            // path: genome fasta file (required)
    gtf              // path: genome GTF file (required)
    additional_fasta // path: additional fasta file to concatenate (optional)

    main:
    ch_versions = Channel.empty()

    //-------------------------------------
    // 1) Process FASTA file (decompress if needed)
    //-------------------------------------
    if (fasta.toString().endsWith('.gz')) {
        ch_fasta_input = Channel.value([[:], fasta])
        GUNZIP_FASTA(ch_fasta_input)
        ch_fasta = GUNZIP_FASTA.out.gunzip.map { meta, file -> file }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(fasta)
    }

    //-------------------------------------
    // 2) Process GTF file (decompress if needed)
    //-------------------------------------
    if (gtf.toString().endsWith('.gz')) {
        ch_gtf_input = Channel.value([[:], gtf])
        GUNZIP_GTF(ch_gtf_input)
        ch_gtf = GUNZIP_GTF.out.gunzip.map { meta, file -> file }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.value(gtf)
    }

    //---------------------------------------------------
    // 3) Handle additional FASTA (if provided)
    //---------------------------------------------------
    ch_final_fasta = ch_fasta
    
    if (additional_fasta) {
        // Process additional FASTA (decompress if needed)
        if (additional_fasta.toString().endsWith('.gz')) {
            ch_additional_input = Channel.value([[:], additional_fasta])
            GUNZIP_ADDITIONAL(ch_additional_input)
            ch_additional_fasta = GUNZIP_ADDITIONAL.out.gunzip.map { meta, file -> file }
            ch_versions = ch_versions.mix(GUNZIP_ADDITIONAL.out.versions)
        } else {
            ch_additional_fasta = Channel.value(additional_fasta)
        }
        
        // For now, we'll use the original fasta
        // TODO: Add concatenation logic when needed
        log.warn "Additional FASTA concatenation not yet implemented. Using original FASTA: ${fasta}"
        // Future implementation: concatenate ch_fasta and ch_additional_fasta
        // ch_final_fasta = concatenate_process(ch_fasta, ch_additional_fasta)
    }

    //----------------------------------------------------
    // 4) Generate STAR index
    //----------------------------------------------------
    STAR_GENOMEGENERATE(
        ch_final_fasta.map { [ [:], it ] },
        ch_gtf.map         { [ [:], it ] }
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    fasta     = ch_final_fasta                               // path: final fasta file
    gtf       = ch_gtf                                       // path: final GTF file  
    index     = STAR_GENOMEGENERATE.out.index.map { it[1] } // path: STAR genome index directory
    versions  = ch_versions                                  // channel: versions.yml files
}
