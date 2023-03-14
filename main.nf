#!/usr/bin/env nextflow

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
    log.info "┇          GeneLab MethylSeq Workflow: $workflow.manifest.version            ┇"
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
if ( ! params.gldsAccession && ! params.runsheet ) {

    println "\n    ${RED}A specific GLDS number (e.g., `--gldsAccession GLDS-397`) or a user-created"
    println "    run sheet (e.g., `--runsheet my-runsheet.csv`) needs to be provided.${NC}"

    println "\n  Exiting for now.\n"
    exit 1

}

/* **** checking that ONLY a glds accession or runsheet was specified **** */
if ( params.gldsAccession && params.runsheet ) {

    println "\n    ${RED}A specific GLDS number (e.g., `--gldsAccession GLDS-397`) AND a user-created"
    println "    run sheet (e.g., `--runsheet my-runsheet.csv`) cannot both be provided.${NC}"

    println "\n  Exiting for now.\n"
    exit 1

}

/* **** making sure if non_directional is set, so is rrbs **** */
// with TRIMGALORE, non_directional can only be specified if it's rrbs (see links in error message below)
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

/* **** checking lib_type set in nextflow.config **** */
if ( params.lib_type !in params.accepted_lib_types ) {

    println "\n    ${RED}No suitable 'lib_type' was set in nextflow.config file.${NC}"

    println "\n  Exiting for now.\n"
    exit 1

}


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
include { TRIMGALORE ; NUGEN_TRIM ; ALIGNMENT_QC } from './modules/QC.nf'
include { PARSE_ANNOTATIONS_TABLE } from './modules/genelab.nf'
include { DOWNLOAD_GUNZIP_REFERENCES ; GTF_TO_PRED ; 
          PRED_TO_BED ; MAKE_GENE_TRANSCRIPT_MAP } from './modules/utilities.nf'
include { GEN_BISMARK_REF ; ALIGN ; DEDUPLICATE ; 
          EXTRACT_METHYLATION_CALLS ; GEN_BISMARK_SAMPLE_REPORT ;
          GEN_BISMARK_SUMMARY } from './modules/bismark.nf'
include { DIFFERENTIAL_METHYLATION_ANALYSIS } from './modules/methylkit.nf'


process TEMP_BLOCK_FOR_PAIRED {

    debug true

    input:
        val(meta)

    output:
        val(meta)

    exec:
        if ( meta.paired_end ) {

            println "\n    ${RED}The workflow is not yet ready to deal with paired-end data, which is what this"
            println "    looks like based on the run sheet.${NC}"

            println "\n  Exiting for now.\n"

            exit 1

        }

}

////////////////////////////////////////////////////
/* --          SUB-WORKFLOWS INCLUDED          -- */
////////////////////////////////////////////////////

include { staging as STAGING } from './staging.nf'


////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////

workflow {

    //staging (and downloading if needed) input reads and metadata
    if ( params.gldsAccession ) { 

        // setting glds_acc channel
        ch_glds_accession = Channel.from( params.gldsAccession )

        // running staging sub-workflow
        STAGING( ch_glds_accession, params.stageLocal )

    } else {

        // setting glds_acc channel to null
        ch_glds_accession = null
        STAGING( ch_glds_accession, params.stageLocal )

    }

    // setting raw reads to channel
    STAGING.out.raw_reads | set { ch_input_reads }

    // storing general info to ch_meta
    STAGING.out.raw_reads | first | 
                            map { it -> it[0] } |
                            // view { meta -> "${YELLOW}  Autodetected Processing Metadata:\n\t pairedEND: ${meta.paired_end}\n\t organism: ${meta.organism_sci}\n\t primary_keytype: ${meta.primary_keytype}${NC}" } |
                            set { ch_meta }

    // adding stop for now if paired-end
    // this PROCESS is defined right above the workflow in this document
    TEMP_BLOCK_FOR_PAIRED( ch_meta ) | set { ch_meta }

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

    // combinging trimgalore logs
    TRIMGALORE.out.reports | map { it -> it[1] } | 
                             collectFile( name: "trimgalore-reports.txt", 
                                          newLine: true, 
                                          storeDir: params.filtered_reads_dir )

    // running NuGEN-specific script if needed
    if ( params.lib_type == 3 ) { 

        ch_nugen_trim_script = channel.fromPath( "bin/trimRRBSdiversityAdaptCustomers.py" )

        NUGEN_TRIM( ch_nugen_trim_script | combine( TRIMGALORE.out.reads ) )

        // setting channel holding nugen-trimmed reads
        ch_trimmed_reads = NUGEN_TRIM.out.reads

        // combining all nugen-trimming logs into one file
        NUGEN_TRIM.out.logs | map { it -> it[1] } |
                              collectFile( name: "nugen-trimming-logs.txt",
                                           newLine: true,
                                           storeDir: params.filtered_reads_dir )

    } else { 

        // updating channel holding initial trimmed reads
            // (so can use same ch variable whether nugen script run or not)
        ch_trimmed_reads = TRIMGALORE.out.reads

    }

    TRIMMED_FASTQC( ch_trimmed_reads )

    // getting all trimmed fastqc output files into one channel
    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } |
                                flatten | collect | set { ch_trimmed_mqc_inputs }

    // multiqc on raw fastqc outputs
    TRIMMED_MULTIQC( ch_trimmed_mqc_inputs )

    // getting target reference info 
    PARSE_ANNOTATIONS_TABLE( params.reference_table_url, ch_meta.organism_sci )

    // downloading reference fasta and gtf files
    DOWNLOAD_GUNZIP_REFERENCES( PARSE_ANNOTATIONS_TABLE.out.reference_genome_urls, PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source )

    // setting reference fasta to channel (for if/when integrating local refs)
    DOWNLOAD_GUNZIP_REFERENCES.out.fasta | set { ch_input_ref }

    // making bismark index from provided input fasta
    GEN_BISMARK_REF( ch_input_ref )

    // aligning 
    ch_trimmed_reads | combine( GEN_BISMARK_REF.out.ch_bismark_index_dir ) | ALIGN

    // combinging aligning reports
    ALIGN.out.reports | map { it -> it[1] } | 
                        collectFile( name: "bismark-align-reports.txt", 
                                     newLine: true, 
                                     storeDir: params.bismark_alignments_dir )

    // Alignment QC
    ALIGNMENT_QC( ALIGN.out.bams, DOWNLOAD_GUNZIP_REFERENCES.out.gtf )

    // deduplicating only if *not* RRBS    
    if ( ! params.rrbs ) {

        DEDUPLICATE( ALIGN.out.bams )

        // setting deduped bams to 'ch_bams_to_extract_from'
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

        // setting channel holding bams to extract from (same as ch_bams_to_extract_from when no dedupe was done)
        ALIGN.out.bams | set { ch_initial_bams }

        // creating empty channel for dedupe reports (so we can pass things the same way to bismark2summary later)
        ALIGN.out.bams | map { it -> [ it[0], '' ] } | set { ch_dedupe_reports }

    }

    // extracting methylation calls
    EXTRACT_METHYLATION_CALLS( ch_bams_to_extract_from, GEN_BISMARK_REF.out.ch_bismark_index_dir )

    // combinging methylation call reports
    EXTRACT_METHYLATION_CALLS.out.reports | map { it -> it[1] } | 
                                            collectFile( name: "bismark-methylation-call-reports.txt", 
                                                         newLine: true, 
                                                         storeDir: params.bismark_methylation_calls_dir )

    // putting all individual sample reports into one channel
    ch_all_sample_reports = ALIGN.out.reports | join( EXTRACT_METHYLATION_CALLS.out.reports ) | 
                                                join( EXTRACT_METHYLATION_CALLS.out.biases ) |
                                                join( ALIGN.out.nuc_stats ) |
                                                join( ch_dedupe_reports )

    // generating individual sample bismark reports
    GEN_BISMARK_SAMPLE_REPORT( ch_all_sample_reports )

    // making channel holding all input files for bismark2summary (bam files, align reports, splitting reports, dedupe reports if they exist)
    ch_bams_and_all_reports = ch_initial_bams | join( ch_all_sample_reports ) | map { it -> it[ 1..it.size() - 1 ] } | collect

    // making overall bismark summary     
    GEN_BISMARK_SUMMARY( ch_bams_and_all_reports )

    // converting GTF to BED
    GTF_TO_PRED( DOWNLOAD_GUNZIP_REFERENCES.out.gtf )
    PRED_TO_BED( GTF_TO_PRED.out.pred )
    
    // making a mapping file of genes to transcripts (needed to link to functional annotations in primary output table)
    MAKE_GENE_TRANSCRIPT_MAP( DOWNLOAD_GUNZIP_REFERENCES.out.gtf )

    // on to R and methylseq next
    // putting runsheet into channel
    ch_runsheet = channel.fromPath( params.runsheet )

    // putting script into channel so works when run from nextflow work dir
    ch_methylkit_script = channel.fromPath( "bin/differential-methylation.R" )

    // putting needed directories in a channel so things can be found when running in the nextflow work dir
    ch_reference_dir = channel.fromPath( params.ref_genome_dir )
    ch_bismark_coverages_dir = channel.fromPath( params.bismark_methylation_calls_dir )

    // need to make coverage files one of the inputs so it knows to wait to start this
    ch_all_bismark_coverage_files = EXTRACT_METHYLATION_CALLS.out.covs | collect
    DIFFERENTIAL_METHYLATION_ANALYSIS( ch_methylkit_script,
                                       ch_all_bismark_coverage_files,
                                       ch_bismark_coverages_dir,
                                       ch_reference_dir,
                                       PARSE_ANNOTATIONS_TABLE.out.simple_organism_name, 
                                       ch_runsheet,
                                       params.reference_table_url,
                                       PARSE_ANNOTATIONS_TABLE.out.annotations_db_url,
                                       ch_meta.primary_keytype )

}
