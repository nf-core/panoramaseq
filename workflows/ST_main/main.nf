/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { FASTQC                 } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// include { paramsSummaryMap       } from 'plugin/nf-schema'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_panoramaseq_pipeline'


include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract'
include { CHECK_SAMPLESHEET } from '../../modules/local/checksamplesheet'
include { CUTADAPT } from '../../modules/nf-core/cutadapt/main'  
include { CUTADAPT as CUTADAPT_ADV } from '../../modules/nf-core/cutadapt/main'
include { FASTQC } from '../../modules/nf-core/fastqc'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'
include { Decode_batch } from '../../modules/local/Decoding/main' 
include { CUTADAPT_ADV_PIPE } from '../../modules/local/cutadapt_adv_pipe/main'
include { STAR_ALIGN_LOCAL } from '../../modules/local/star_align/main'
include { SAMTOOLS_INDEX as index1; SAMTOOLS_INDEX as index2 } from '../../modules/nf-core/samtools/index/main'
include { FEATURECOUNTS_CUSTOM } from '../../modules/local/FeatCounts/main'   
include { UMI_count } from '../../modules/local/UMI_count/main'
include { SAMTOOLS_SORT_LOCAL } from '../../modules/local/SAMTOOLS_SORT_LOCAL/main' 
// include { MULTIQC } from '../../modules/nf-core/multiqc/main'
// include { PIPELINE_COMPLETION } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ST_main {
    
    take:

    valid_data // channel: samplesheet read in from checked samplesheet process
    
    main:

    ch_versions = Channel.empty()
    
    // 2. Subsample reads for samples with sample_size set in meta using SEQTK_SAMPLE
    seqtk_input = valid_data.filter{ meta, reads -> meta.sample_size != null }
    SEQTK_SAMPLE(seqtk_input)

    // 3. For samples without sample_size, pass through directly
    passthrough = valid_data.filter{ meta, reads -> meta.sample_size == null }

    // 4. Combine both channels for downstream processing
    // Use concat instead of mix to avoid multi-channel operator issues
    umitools_input = SEQTK_SAMPLE.out.reads.concat(passthrough)

    // 4. Extract UMIs using UMITOOLS_EXTRACT
    umi_extract = UMITOOLS_EXTRACT(umitools_input)

    // 5. Trim reads after UMI extraction using CUTADAPT
    cutadapt_results = CUTADAPT(UMITOOLS_EXTRACT.out.reads)

    // 6. Decode barcodes using Decode_batch
    decode_results = Decode_batch(CUTADAPT.out.reads)

    // 7. Advanced trimming on R2 using CUTADAPT_ADV_PIPE
    cutadapt2_results = CUTADAPT_ADV_PIPE(Decode_batch.out.reads)

    // 8. Align single-end trimmed R2 fastq using STAR_ALIGN_LOCAL
    //     Requires reference genome directory
    //     Exits if not specified
    //     Prepares input as tuple of meta, fastq, and genome dir
    def star_genome_dir = params.star_genome_dir
    if (!star_genome_dir) {
        exit 1, "ERROR: STAR reference genome directory must be specified via --star_genome_dir."
    }
    star_local_input = CUTADAPT_ADV_PIPE.out.reads.map { meta, r2fastq ->
        [meta, file(r2fastq), file(star_genome_dir)]
    }
    STAR_ALIGN_LOCAL(star_local_input)

    // 9. Index the sorted BAM output from STAR_ALIGN_LOCAL using samtools index (index1)
    samtools_index_input = STAR_ALIGN_LOCAL.out.bam
    samtools_index_results = samtools_index_input | index1

    // 10. Count features from BAM using annotation with FEATURECOUNTS_CUSTOM
    //     Exits if annotation file is not specified
    //     Prepares input as tuple of meta, bam, and annotation file
    def annotation_file = params.star_gtf
    if (!annotation_file) {
        exit 1, "ERROR: GTF annotation file must be specified via --star_gtf."
    }
    custom_featurecounts_input = STAR_ALIGN_LOCAL.out.bam.map { meta, bam ->
        tuple(meta, bam, file(annotation_file))
    }
    FEATURECOUNTS_CUSTOM(custom_featurecounts_input)

    // 11. Sort BAM files after feature counting using SAMTOOLS_SORT_LOCAL
    SAMTOOLS_SORT_LOCAL(FEATURECOUNTS_CUSTOM.out.annotated_bam.map { meta, bam -> tuple(meta, bam) })

    // 12. Index the sorted BAM output from SAMTOOLS_SORT_LOCAL using samtools index (index2)
    samtools_index2_input = SAMTOOLS_SORT_LOCAL.out.bam.map { meta, bam -> tuple(meta, bam) }
    samtools_index2_results = samtools_index2_input | index2

    // 13. Count UMIs using UMI_count
    //     Joins sorted BAM and its index, then passes as tuple to UMI_count
    UMI_count_input = SAMTOOLS_SORT_LOCAL.out.bam.join(index2.out.bai)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }
    UMI_count(UMI_count_input)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
