/*
 * Processes related to bismark processing.
 */

process GEN_BISMARK_REF {

    debug true
    // tag "On: $name"

    // input:
    //     tuple val(name), path(reads)

    // output:
    //     tuple val(name), path("${ name }${ params.file_suffix }_fastqc.zip"), path("${name}${ params.file_suffix }_fastqc.html"), emit: fastqc

    script:

        """
        bismark -v
        bowtie2 --version
        samtools --version        
        """

}
