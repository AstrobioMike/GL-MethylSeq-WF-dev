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


// process basicExample {
    
//     input:
//         val x from num
    
//     "echo process job $x"

// }

////////////////////////////////////////////////////
/* --            STARTING CHANNELS             -- */
////////////////////////////////////////////////////

ch_x = channel.fromList([1, 2])
ch_y = channel.fromList(['a', 'b', 'c'])


////////////////////////////////////////////////////
/* --                PROCESSES                 -- */
////////////////////////////////////////////////////


process basicExample {
    
    input:
        val num

    output:
        stdout

    "echo process job $num"

}


process REPORT {

    name = "REPORT"

    input:
        path input_reads

    output:
        stdout emit: stdout

    script:
        """
        echo "Doing stuff to $input_reads"
        echo "Process task name is $task.name"
        ls -l $input_reads
        """
}


process tupleExample {

    input:
        tuple val(x), path('latin.txt')

    """
    echo Processing $x
    cat - latin.txt > copy
    """

}


process foo {

    debug true

    input:
        val x
        val y

    script:
        """
        echo $x and $y
        """

}

workflow {

    // num = channel.from( 1, 2, 3 )
    // input_reads = channel.fromPath(params.input_reads)
    // values = channel.from( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

    // // REPORT(input_reads)
    // // REPORT.out.stdout | view

    // // basicExample(num) | view

    // tupleExample(values)

    foo( ch_x, ch_y )

}
