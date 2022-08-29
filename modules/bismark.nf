/*
 * Processes related to bismark processing.
 */

process GEN_BISMARK_REF {


    // I'm trying to get all the outputs into ./Bismark_Index
    // these should include: "bismark-.*.txt", "tiny-ref-fasta.fa", and "Bisulfite_Genome/*"
    
    // if i just have this, it, understandably, puts all the outputs except the "bismark-.*.txt" stdout/err capture file there:
    publishDir "./", mode: 'link'

    // but once i tried also sending the stdout/err there i couldn't figure out how to recover/specify all the other outputs anymore
    // e.g. this doesn't work:
    // publishDir "./", mode: 'link', pattern: "${ params.bismark_index_dir }**"
    // publishDir params.bismark_index_dir, mode: 'link', pattern: "bismark-genome-preparation-output.txt"
    // as the only thing that ends up in the right place is the stdout/err file
    // so i think my pattern isn't working properly

    // Do you spot what i can do? If this makes any sense...
    // It might look a little extra odd because of how bismark_genome_preparation wants to act on a directory, and how i've integrated that here

    input:
        path( input_ref_fasta )

    output:
        path( params.bismark_index_dir ), emit: ch_bismark_index_dir
        path("bismark-genome-preparation-output.txt")

    script:

        """
        mkdir -p ${params.bismark_index_dir}

        cp ${input_ref_fasta} ${params.bismark_index_dir}

        bismark_genome_preparation --bowtie2 --parallel ${params.bismark_index_creation_threads} ${params.bismark_index_dir} > bismark-genome-preparation-output.txt 2>&1
        """

}
