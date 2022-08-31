/*
 * Processes related to bismark processing.
 */

process GEN_BISMARK_REF {

    publishDir "./", mode: 'link'

    input:
        path( input_ref_fasta )

    output:
        path( params.bismark_index_dir ), emit: ch_bismark_index_dir

    script:

        """
        mkdir -p ${params.bismark_index_dir}

        # copying ref fasta in there so everything exists with the bismark index dir (required by bismark)
        cp ${input_ref_fasta} ${params.bismark_index_dir}

        # making index
        bismark_genome_preparation --bowtie2 --parallel ${params.bismark_index_creation_threads} ${params.bismark_index_dir} > bismark-genome-preparation-output.txt 2>&1

        # moving log file into bismark dir since i can't figure out how to do it with publishDirs above
        mv bismark-genome-preparation-output.txt ${params.bismark_index_dir}
        """        

}


// process ALIGN {

//     publishDir params.bismark_alignments, mode: 'link'

//     input:
//         path( reads )
//         path( bismark_index_dir )

//     output:
//         path( params.bismark_index_dir ), emit: ch_bismark_index_dir
//         // path("bismark-genome-preparation-output.txt")

//     script:

//         """
//         mkdir -p ${params.bismark_index_dir}

//         # copying ref fasta in there so everything exists with the bismark index dir (required by bismark)
//         cp ${input_ref_fasta} ${params.bismark_index_dir}

//         # making index
//         bismark_genome_preparation --bowtie2 --parallel ${params.bismark_index_creation_threads} ${params.bismark_index_dir} > bismark-genome-preparation-output.txt 2>&1

//         # moving log file into bismark dir since i can't figure out how to do it with publishDirs above
//         mv bismark-genome-preparation-output.txt ${params.bismark_index_dir}
//         """        

// }
