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


process ALIGN {

    tag "On: $name"

    publishDir params.bismark_alignments_dir, mode: 'link', pattern: "${ name }_trimmed_bismark_*.bam"

    input:
        tuple val(name), path(reads), path( bismark_index_dir )

    output:
        tuple val(name), path("${ name }_trimmed_bismark_*.bam"), emit: bams
        path("${ name }_trimmed_bismark_*_report.txt"), emit: reports

    script:

        non_directional = params.non_directional ? '--non_directional' : ''
        fastq_files = params.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

        """
        bismark --bam -p ${params.bismark_align_threads} --genome ${bismark_index_dir} ${fastq_files}
        """

}


// process DEDUPLICATE {



// }