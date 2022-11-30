#!/usr/bin/env nextflow
/*
==========================================================================================
Largely modified from the nf-core/methylseq workflow: https://github.com/nf-core/methylseq
==========================================================================================
*/

// Declare syntax version
nextflow.enable.dsl=2

// making sure specific expected nextflow version is being used
if( ! nextflow.version.matches( workflow.manifest.nextflowVersion ) ) {
    println "\n    ${RED}This workflow requires Nextflow version $workflow.manifest.nextflowVersion, but version $nextflow.version is currently active.${NC}"
    println "\n    ${YELLOW}You can set the proper version for this terminal session by running: `export NXF_VER=$workflow.manifest.nextflowVersion`${NC}"

    println "\n  Exiting for now.\n"
    exit 1
}


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

if (params.help) {
    log.info "\n┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅"
    log.info "┇          GeneLab Methyl-seq Workflow: $workflow.manifest.version            ┇"
    log.info "┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅\n"
    
    log.info "                        HELP MENU\n"

    log.info "General settings can be configured in the 'nextflow.config' file.\n"

    log.info "  ${YELLOW}Usage example 1: Processing a GLDS dataset${NC}\n"
    log.info "    `nextflow run main.nf --gldsAccession GLDS-397`\n"

    log.info "  ${YELLOW}Usage example 2: Processing a local dataset (requires user-created runsheet)${NC}\n"
    log.info "    `nextflow run main.nf --runsheet my-runsheet.csv`\n\n"

    exit 0
}


////////////////////////////////////////////////////
/* --            PRE-FLIGHT CHECKS             -- */
////////////////////////////////////////////////////

/* **** checking a glds accession or runsheet was specified **** */
if ( ! params.gldsAccession ) {

    if ( ! params.runsheet ) {

        println "\n    ${RED}A specific GLDS number (e.g., `--gldsAccession GLDS-397`) or a user-created"
        println "    run sheet (e.g., `--runsheet my-runsheet.csv`) needs to be provided.${NC}"

        println "\n  Exiting for now.\n"
        exit 1

    }

}

/* **** making sure if non_directional is set, so is rrbs **** */
// with TRIMGALORE, non_directional can only be specified if it's rrbs
// https://github.com/FelixKrueger/TrimGalore/blob/e4fb81ff4adf052cd22e859f0f36aaee7ce63489/trim_galore#L2637
// https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#non-directional-mode
// see last 2 slides here: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf
// so, making sure that is the case here (might have to adjust if we auto-set the non_directional flag)
if ( params.non_directional == true ) {

    if ( params.rrbs != true ) {

        println "\n    ${RED}As per TrimGalore, the non_directional parameter can only be set to 'true' if the rrbs parameter is also set to 'true'.${NC}"
        println "\n    ${YELLOW}See:"
        println "        - https://github.com/FelixKrueger/TrimGalore/blob/e4fb81ff4adf052cd22e859f0f36aaee7ce63489/trim_galore#L2637"
        println "        - https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#non-directional-mode"
        println "        - and last few slides here: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf${NC}"

        println "\n  Exiting for now.\n"
        exit 1

    }
}




/* **** checking specified target_organism is available **** */
// if ( params.target_organism !in params.accepted_target_orgs ) {

//     println "\n    ${RED}No suitable 'target_organism' was set in nextflow.config file.${NC}"
//     println "\n    Currently available options are:\n"
//     for ( org in params.accepted_target_orgs ) {

//         println "        ${org}"

//     }

//     println "\n  Exiting for now.\n"

//     exit 1

// }

/* **** checking lib_type set in nextflow.config **** */
if ( params.lib_type !in params.accepted_lib_types ) {

    println "\n    ${RED}No suitable 'lib_type' was set in nextflow.config file.${NC}"

    println "\n  Exiting for now.\n"
    exit 1

}

/* **** checking specified input_reads_dir exists **** */
// if ( ! file( params.input_reads_dir ).exists() ) {

//     println "\n    ${RED}The specified '${params.input_reads_dir}' directory set in nextflow.config can't be found.${NC}"
//     println "\n  Exiting for now.\n"

//     exit 1

// }

/* **** checking for gzipped reads (anything ending with fq.gz or fastq.gz) **** */
// creating and adding files to list
// input_file_list = []
// file(params.input_reads_dir).eachFileMatch(~/.*fastq.gz|.*.fq.gz/) { target_file ->

//     input_file_list << target_file

// }

// exiting and reporting if none found
// if ( input_file_list.size() == 0 ) {

//     println "\n    ${RED}No gzipped fastq files were found in the specified ${params.input_reads_dir} directory set in nextflow.config.${NC}"
//     println "  Exiting for now.\n"

//     exit 1

// }

// Right now i can't get fastqc/multiqc to retain things like "_trimmed" in the filenames, so
// things are getting collapsed if the input files don't end with "_raw".
// Putting in a stop-gap for now that exits if the input files don't include "_raw":
// for ( input_file in input_file_list ) {

//     if ( ! input_file.toString().endsWith("_raw.fastq.gz") ) {
        
//         if ( ! input_file.toString().endsWith("_raw.fq.gz") ) {

//             println "\n    ${RED}Currently the input read files need to have a suffix like '_raw.fastq.gz' or '_raw.fq.gz'.${NC}"
//             println "  Exiting for now.\n"

//             exit 1

//         }

//     }

// }


////////////////////////////////////////////////////
/* --              DEBUG WARNING               -- */
////////////////////////////////////////////////////
if ( params.truncateTo || ! params.stageLocal ) {

    println( "\n${YELLOW}    WARNING WARNING: DEBUG OPTIONS ARE ENABLED!${NC}\n" )

    params.limitSamplesTo ? println( "${YELLOW}        - Limiting run to ${ params.limitSamplesTo } samples${NC}" ) : null
    params.truncateTo ? println( "${YELLOW}        - Truncating reads to first ${ params.truncateTo } records${NC}" ) : null
    params.genomeSubsample ? println( "${YELLOW}        - Subsetting reference genome to just the '${ params.genomeSubsample }' contig${NC}" ) : null
    params.stageLocal ? null : println( "${YELLOW}        - stageLocal has been set to false (reads won't be downloaded)${NC}" )

    println ""

}


////////////////////////////////////////////////////
/* --            PROCESSES INCLUDED            -- */
////////////////////////////////////////////////////

include { FASTQC as RAW_FASTQC } from './modules/QC.nf' addParams( file_suffix: "" )
include { FASTQC as TRIMMED_FASTQC } from './modules/QC.nf' addParams( file_suffix: "_trimmed" )
include { MULTIQC as RAW_MULTIQC } from './modules/QC.nf' addParams( MQCLabel: "raw" )
include { MULTIQC as TRIMMED_MULTIQC } from './modules/QC.nf' addParams( MQCLabel: "trimmed" )
include { MULTIQC as PROJECT_MULTIQC } from './modules/QC.nf' addParams( MQCLabel: "project" )
include { TRIMGALORE ; ALIGNMENT_QC } from './modules/QC.nf'
include { PARSE_ANNOTATIONS_TABLE } from './modules/genelab.nf'
include { DOWNLOAD_GUNZIP_REFERENCES ; GTF_TO_PRED ; 
          PRED_TO_BED ; MAKE_GENE_TRANSCRIPT_MAP } from './modules/utilities.nf'
include { GEN_BISMARK_REF ; ALIGN ; DEDUPLICATE ; 
          EXTRACT_METHYLATION_CALLS ; GEN_BISMARK_SAMPLE_REPORT ;
          GEN_BISMARK_SUMMARY } from './modules/bismark.nf'


////////////////////////////////////////////////////
/* --          SUB-WORKFLOWS INCLUDED          -- */
////////////////////////////////////////////////////

include { staging as STAGING } from './staging.nf'
// include { references as REFERENCES } from './references.nf'


////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////

workflow {

    /* **** staging (and downloading if needed) input reads and metadata **** */

    if ( params.gldsAccession ) { 

        // setting glds_acc channel
        ch_glds_accession = Channel.from( params.gldsAccession)

        // running staging sub-workflow
        STAGING( ch_glds_accession, params.stageLocal )

    } else {

        println "\n    ${RED}We're not ready to take a user-provided run sheet yet :(${NC}"
        println "\n  Exiting for now.\n"

        exit 1

    }

    // setting raw reads to channel
    STAGING.out.raw_reads | set { ch_input_reads }

    // storing general info to ch_meta
    STAGING.out.raw_reads | first | 
                            map { it -> it[0] } |
                            view { meta -> "${YELLOW}  Autodetected Processing Metadata:\n\t pairedEND: ${meta.paired_end}\n\t organism: ${meta.organism_sci}${NC}" } |
                            set { ch_meta }

    // // // detecting input reads and removing extensions from their unique sample names
    // // ch_input_reads = Channel.fromFilePairs( input_file_list, size: params.single_end ? 1 : 2 ) { file -> file.name.replaceAll( /.fastq.gz|.fq.gz|_raw.fastq.gz|_raw.fq.gz/, '' ) }
    
    // // // writing out unique sample names to file and setting to channel
    // // ch_input_reads | map { it -> it[0] } |
    // //                  collectFile( name: 'samples.txt', newLine: true, storeDir: "./" ) |
    // //                  set { ch_samples_txt }

    

    // raw fastqc on input reads
    RAW_FASTQC( ch_input_reads )
    
    // getting all raw fastqc output files into one channel
    RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2]] } |
                            flatten | collect | set { ch_raw_mqc_inputs }

    // setting channel for multiqc config file
    ch_multiqc_config = Channel.fromPath( params.multiqc_config )

    // multiqc on raw fastqc outputs
    RAW_MULTIQC( ch_raw_mqc_inputs )

    // quality trimming/filtering input reads
    TRIMGALORE( ch_input_reads )

    // combinging trimming logs
    TRIMGALORE.out.reports | map { it -> it[1] } | 
                             collectFile( name: "trimgalore-reports.txt", 
                                          newLine: true, 
                                          storeDir: params.filtered_reads_dir )

    // // fastqc on trimmed reads
    TRIMMED_FASTQC( TRIMGALORE.out.reads )

    // getting all trimmed fastqc output files into one channel
    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } |
                                flatten | collect | set { ch_trimmed_mqc_inputs }

    // multiqc on raw fastqc outputs
    TRIMMED_MULTIQC( ch_trimmed_mqc_inputs )

    // // // setting input reference fasta file channel
    // // ch_input_ref = Channel.fromPath( params.genome, checkIfExists: true )

    // getting target reference info 
    PARSE_ANNOTATIONS_TABLE( params.reference_table_url, ch_meta.organism_sci )

    // downloading reference fasta and gtf files
    DOWNLOAD_GUNZIP_REFERENCES( PARSE_ANNOTATIONS_TABLE.out.reference_genome_urls, PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source )

    // setting reference fasta to channel (for if/when integrating local refs)
    DOWNLOAD_GUNZIP_REFERENCES.out.fasta | set { ch_input_ref }

    // making bismark index from provided input fasta
    GEN_BISMARK_REF( ch_input_ref )

    // aligning 
    TRIMGALORE.out.reads | combine( GEN_BISMARK_REF.out.ch_bismark_index_dir ) | ALIGN

    // combinging aligning reports
    ALIGN.out.reports | map { it -> it[1] } | 
                        collectFile( name: "bismark-align-reports.txt", 
                                     newLine: true, 
                                     storeDir: params.bismark_alignments_dir )

    // deduplicating only if *not* RRBS    
    if ( ! params.rrbs ) {

        DEDUPLICATE( ALIGN.out.bams )

        // setting deduped bams to 'ch_bams'
        DEDUPLICATE.out.bams | set { ch_bams_to_extract_from }

        // setting channel holding initial bams
        ALIGN.out.bams | set { ch_initial_bams }

        // combining dedupe reports
        DEDUPLICATE.out.reports | map { it -> it[1] } | 
                                  collectFile( name: "bismark-dedupe-reports.txt", 
                                               newLine: true, 
                                               storeDir: params.bismark_alignments_dir )

        // setting deduped reports channel to 'ch_dedupe_reports'
        DEDUPLICATE.out.reports | set { ch_dedupe_reports }
    

    } else {

        // setting channel holding initial bams
        ALIGN.out.bams | set { ch_bams_to_extract_from }

        // creating empty channel for dedupe reports and deduped bams
        ALIGN.out.bams | map { it -> [ it[0], '' ] } | set { ch_dedupe_reports }

    }

    // extracting methylation calls
    EXTRACT_METHYLATION_CALLS( ch_bams_to_extract_from )

    // combinging methylation call reports
    EXTRACT_METHYLATION_CALLS.out.reports | map { it -> it[1] } | 
                                            collectFile( name: "bismark-methylation-call-reports.txt", 
                                                         newLine: true, 
                                                         storeDir: params.bismark_methylation_calls_dir )

    // putting all individual sample reports into one channel
    ch_all_sample_reports = ALIGN.out.reports | join( EXTRACT_METHYLATION_CALLS.out.reports ) | 
                                                join( EXTRACT_METHYLATION_CALLS.out.biases ) |
                                                join( ch_dedupe_reports )

    // generating individual sample bismark reports
    GEN_BISMARK_SAMPLE_REPORT( ch_all_sample_reports )

    // making channel holding all input files for bismark2summary (bam files, align reports, splitting reports, dedupe reports)
    ch_bams_and_all_reports = ch_bams_to_extract_from | join( ch_all_sample_reports ) | map { it -> it[ 1..it.size() - 1 ] } | collect

    // making overall bismark summary 
        // problem with this for now, see issue i posted here: https://github.com/FelixKrueger/Bismark/issues/520
            // ahh, bismark2summary needs the original bams to start with even when deduplicated (makes sense), passing them too now
    // so if rrbs, then adding initial bams to ch_bams_and_all_reports
    if ( ! params.rrbs ) {

        // in here, the map{} is to drop the meta tags so that the mix works properly (i think needed?)
        ch_initial_bams = ch_initial_bams | map { it -> it[1] } | collect

        ch_bams_and_all_reports = ch_bams_and_all_reports | mix ( ch_initial_bams ) | collect

    } 
    
    GEN_BISMARK_SUMMARY( ch_bams_and_all_reports )

    // Alignment QC
    ALIGNMENT_QC( ch_bams_to_extract_from )

    // generate multiqc project report
        // creating input channel holding all needed inputs for the project-level multiqc

        // passing the projectDir variable as a channel to grab everything
    full_project_dir_ch = Channel.fromPath( projectDir )

        // adding additional needed channels
    full_project_dir_ch | mix( ALIGN.out.reports |  map { it -> it[1] }, 
                               EXTRACT_METHYLATION_CALLS.out.reports | map { it -> it[1] },
                               ch_raw_mqc_inputs,
                               ch_trimmed_mqc_inputs,
                               ALIGNMENT_QC.out.qualimaps
                             ) | 
                          collect | set{ project_multiqc_in_ch }


    PROJECT_MULTIQC( project_multiqc_in_ch )

    // converting GTF to BED
    GTF_TO_PRED( DOWNLOAD_GUNZIP_REFERENCES.out.gtf )
    PRED_TO_BED( GTF_TO_PRED.out.pred )
    
    // making a mapping file of genes to transcripts (needed to link to functional annotations in primary output table)
    MAKE_GENE_TRANSCRIPT_MAP( DOWNLOAD_GUNZIP_REFERENCES.out.gtf )

    // on to R next... 
    // here is how to pass the link to R script: 
        // PARSE_ANNOTATIONS_TABLE.out.annotations_db_url
    // look at ch_meta for primary_keytype to try to find keytype to pass here
    // figure out what we want to pass for the "assembly" slot in methylkit object,
        // maybe something from ensemblVersion/Source in PARSE_ANNOTATIONS_TABLE.out also


}
