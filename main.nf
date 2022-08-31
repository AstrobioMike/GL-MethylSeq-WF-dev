#!/usr/bin/env nextflow
/*
==========================================================================================
Largely modified from the nf-core/methylseq workflow: https://github.com/nf-core/methylseq
==========================================================================================
*/

// Declare syntax version
nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

if (params.help) {
    log.info "\n┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅"
    log.info "┇          GeneLab Methyl-seq Workflow: $workflow.manifest.version            ┇"
    log.info "┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅\n"
    log.info "    Usage example (after setting parameters in the 'nextflow.config' file):"
    log.info "        `nextflow run main.nf`\n"
    exit 0
}

////////////////////////////////////////////////////
/* --            PRE-FLIGHT CHECKS             -- */
////////////////////////////////////////////////////

/* **** checking lib_type set in nextflow.config **** */
if ( params.lib_type !in params.accepted_lib_types ) {

    println "\n    ${RED}No suitable 'lib_type' was set in nextflow.config.${NC}"
    println "    Exiting for now.\n"

    exit 1

}

/* **** checking specified input_reads_dir exists **** */
if ( ! file( params.input_reads_dir ).exists() ) {

    println "\n    ${RED}The specified '${params.input_reads_dir}' directory set in nextflow.config can't be found.${NC}"
    println "    Exiting for now.\n"

    exit 1

}

/* **** checking for gzipped reads (anything ending with fq.gz or fastq.gz) **** */
// creating and adding files to list
input_file_list = []
file(params.input_reads_dir).eachFileMatch(~/.*fastq.gz|.*.fq.gz/) { target_file ->

    input_file_list << target_file

}

// exiting and reporting if none found
if ( input_file_list.size() == 0 ) {

    println "\n    ${RED}No gzipped fastq files were found in the specified ${params.input_reads_dir} directory set in nextflow.config.${NC}"
    println "    Exiting for now.\n"

    exit 1

}


////////////////////////////////////////////////////
/* --            PROCESSES INCLUDED            -- */
////////////////////////////////////////////////////

include { FASTQC as RAW_FASTQC } from './modules/QC.nf' addParams( file_suffix: "" )
include { FASTQC as TRIMMED_FASTQC } from './modules/QC.nf' addParams( file_suffix: "_trimmed" )
include { MULTIQC as RAW_MULTIQC } from './modules/QC.nf' addParams( MQCLabel: "raw" )
include { MULTIQC as TRIMMED_MULTIQC } from './modules/QC.nf' addParams( MQCLabel: "trimmed" )
include { TRIMGALORE } from './modules/QC.nf'
include { GEN_BISMARK_REF ; ALIGN ; DEDUPLICATE ; EXTRACT_METHYLATION_CALLS } from './modules/bismark.nf'

////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////


workflow {

    // detecting input reads and removing extensions from their unique sample names
    ch_input_reads = Channel.fromFilePairs( input_file_list, size: params.single_end ? 1 : 2 ) { file -> file.name.replaceAll( /.fastq.gz|.fq.gz/,'' ) }

    // writing out unique sample names to file and setting to channel
    ch_input_reads | map { it -> it[0] } |
                     collectFile( name: 'samples.txt', newLine: true, storeDir: "./" ) |
                     set { ch_samples_txt }

    // raw fastqc on input reads
    RAW_FASTQC( ch_input_reads )

    // getting all raw fastqc output files into one channel
    RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2]] } |
                            flatten | collect | set { ch_raw_mqc_inputs }

    // multiqc on raw fastqc outputs
    RAW_MULTIQC( ch_raw_mqc_inputs )

    // quality trimming/filtering input reads
    TRIMGALORE( ch_input_reads )

    // combinging trimming logs
    TRIMGALORE.out.reports | collectFile( name: "trimgalore-reports.txt", 
                                          newLine: true, 
                                          storeDir: params.filtered_reads_dir )

    // fastqc on trimmed reads
    TRIMMED_FASTQC( TRIMGALORE.out.reads )

    // getting all trimmed fastqc output files into one channel
    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2]] } |
                                flatten | collect | set { ch_trimmed_mqc_inputs }

    // multiqc on raw fastqc outputs
    TRIMMED_MULTIQC( ch_trimmed_mqc_inputs )

    // setting input reference fasta file channel
    ch_input_ref = Channel.fromPath( params.genome, checkIfExists: true )

    // making bismark index    
    GEN_BISMARK_REF( ch_input_ref )

    // aligning 
    TRIMGALORE.out.reads | combine( GEN_BISMARK_REF.out.ch_bismark_index_dir ) | ALIGN

    // combinging aligning reports
    ALIGN.out.reports | collectFile( name: "bismark-align-reports.txt", 
                                     newLine: true, 
                                     storeDir: params.bismark_alignments_dir )

    // deduplicating only if *not RRBS    
    if ( ! params.rrbs ) {

        DEDUPLICATE( ALIGN.out.bams )

        // setting deduped bams to 'ch_bams'
        DEDUPLICATE.out.bams | set { ch_bams }

    } else {

        // setting non-deduped bams to 'ch_bams'
        ALIGN.out.bams | set { ch_bams }

    }

    // extracting methylation calls
    EXTRACT_METHYLATION_CALLS( ch_bams )

    // combinging methylation call reports
    EXTRACT_METHYLATION_CALLS.out.reports | collectFile( name: "bismark-methylation-call-reports.txt", 
                                                         newLine: true, 
                                                         storeDir: params.bismark_methylation_calls_dir)

}
