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
/* --                PROCESSES                 -- */
////////////////////////////////////////////////////


// process basicExample {
    
//     input:
//         val num

//     output:
//         stdout

//     "echo process job $num"

// }


// process REPORT {

//     name = "REPORT"

//     input:
//         path input_reads

//     output:
//         stdout emit: stdout

//     script:
//         """
//         echo "Doing stuff to $input_reads"
//         echo "Process task name is $task.name"
//         ls -l $input_reads
//         """
// }


// process tupleExample {

//     input:
//         tuple val(x), path('latin.txt')

//     """
//     echo Processing $x
//     cat - latin.txt > copy
//     """

// }


// process foo {

//     debug true

//     input:
//         val x
//         val y

//     script:
//         """
//         echo $x and $y
//         """

// }

process fastqc {

    tag "$name"

    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"

    input:
        tuple val(name), path(reads)

    output:
        path "*_fastqc.{zip,html}"
        

    script:

        """
        fastqc $reads
        """

}

workflow {

    input_reads = Channel.fromFilePairs( params.input_reads, size: params.single_end ? 1 : 2 ) { file -> file.name.replaceAll( /.fastq.gz|.fq.gz/,'' ) }
    // input_reads.view()

    fastqc(input_reads)

    // num = channel.from( 1, 2, 3 )
    // input_reads = channel.fromPath(params.input_reads)
    // values = channel.from( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

    // // REPORT(input_reads)
    // // REPORT.out.stdout | view

    // // basicExample(num) | view

    // tupleExample(values)

    // foo( ch_x, ch_y )

}
