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