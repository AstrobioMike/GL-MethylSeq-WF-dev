/*
* Workflow taken/modified from Jonathan Oribello: https://github.com/asaravia-butler/GeneLab_Data_Processing/blob/b6576eade4b77e3bbf5bcad573c8d77fa7d20d1c/RNAseq/Workflow_Documentation/NF_RCP-F/workflow_code/NF_RCP-F_1.0.0/stage_analysis.nf
*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates RNASeq Samplesheet
* 2a. Downloads Raw Reads
* 2b. Downloads Truncated Raw Reads (Useful for testing with limited resources)
*/

// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --            PROCESSES INCLUDED            -- */
////////////////////////////////////////////////////

include { METHYLSEQ_RUNSHEET_FROM_GLDS as GENERATE_RUNSHEET;
          GENERATE_METASHEET;
          STAGE_RAW_READS; get_runsheet_paths } from './modules/genelab.nf'


////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////

workflow staging {

    take:
        ch_glds_accession
        stageLocal
    
    main:
        sample_limit = params.limitSamplesTo ? params.limitSamplesTo : -1 // -1 in take means no limit

        if ( ! params.runsheet ) {

            ch_glds_accession | GENERATE_RUNSHEET
        
            GENERATE_RUNSHEET.out.runsheet | set{ ch_runsheet }
        
            GENERATE_METASHEET( GENERATE_RUNSHEET.out.isazip, GENERATE_RUNSHEET.out.runsheet )

        } else {

            ch_runsheet = channel.fromPath(params.runsheet)

        }


        ch_runsheet | splitCsv(header: true)
                    | map{ row -> get_runsheet_paths(row) }
                    | map{ it -> params.force_single_end ? mutate_to_single_end(it) : it }
                    | take( sample_limit )
                    | set{ ch_samples }


        if ( stageLocal && params.truncateTo ) {

            // download truncated raw reads
            // download full raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]] }
                       | branch {
                            paired: it.size() == 3
                            single: it.size() == 2
                         }
                       | set{ ch_raw_read_pointers }


            // PAIRED END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                        | map { it -> [ it[0], [ it[1], it[2] ] ] }
                                        // | view { it -> "TRUNCATED PAIRED READS ($params.truncateTo): $it[0]"}
                                        | set { ch_raw_reads }
            // SINGLE END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                        | map { it -> [ it[0], [ it[1] ] ] }
                                        // | view { it -> "TRUNCATED SINGLE READS ($params.truncateTo): $it[0]"}
                                        | mix( ch_raw_reads )
                                        | set { ch_raw_reads }

            // Moves the truncated files to expected raw read locations as per samplesheet
            ch_raw_reads | STAGE_RAW_READS

        } else if ( stageLocal && ! params.truncateTo ) {

            // download full raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                       | set { ch_raw_reads }

            // Download the raw reads and publish them to expected raw read locations as per samplesheet
            ch_raw_reads | STAGE_RAW_READS

        } else {
      
            // Don't download any raw reads
        }

    emit:
        raw_reads = stageLocal ? STAGE_RAW_READS.out : null
        isa = params.runsheet ? null : GENERATE_RUNSHEET.out.isazip
        runsheet = ch_runsheet
        metasheet = params.runsheet ? null : GENERATE_METASHEET.out.metasheet

}
