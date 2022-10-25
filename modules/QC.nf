/*
 * Processes related to quality assessment/quality control.
 */

process FASTQC {

    tag "On: $meta.id"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${ meta.id }*${ params.file_suffix }_fastqc.zip"), path("${ meta.id }*${ params.file_suffix }_fastqc.html"), emit: fastqc

    script:

        """
        fastqc $reads
        """

}


process MULTIQC {

    tag "On: ${ params.MQCLabel }"

    publishDir params.multiqc_outputs_dir, mode: 'link'

    input:
        // path("samples.txt")
        path("mqc_in/*") // any number of multiqc compatible files
        // path(multiqc_config)

    output:
        path("${ params.MQCLabel }_multiqc.html"), emit: html
        path("${ params.MQCLabel }_multiqc_report.zip"), emit: zipped_report

    script:

        """
        multiqc --interactive -o ${ params.MQCLabel }_multiqc_report \
                -n ${ params.MQCLabel }_multiqc \
                --config ${ params.multiqc_config } mqc_in/*

        mv ${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html .
        
        zip -m -r '${ params.MQCLabel }_multiqc_report.zip' '${ params.MQCLabel }_multiqc_report'
        """

}


process TRIMGALORE {

    tag "On: $meta.id"

    publishDir params.filtered_reads_dir, mode: 'link', pattern: "${ meta.id }_trimmed.*.gz"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${ meta.id }_trimmed.*.gz"), emit: reads
        tuple val(meta), path("${ meta.id }*trimming_report.txt"), emit: reports

    script:

        non_directional = params.non_directional ? '--non_directional' : ''

        """
        # this depends on the lib_type and then if paired-end or not
        if [ ${params.lib_type} == 1 ]; then

            if [ ${params.single_end} == 'true' ]; then

                # trimming
                trim_galore --cores 4 --gzip $reads ${non_directional}

                # renaming to our convention
                mv ${meta.id}*_trimmed.fq.gz ${meta.id}_trimmed.fastq.gz

            else

                printf "    first lib type, paired-end not yet setup!\n"
                exit 1

            fi
        
        elif [ ${params.lib_type} == 2 ]; then

            printf "    second lib type\n"

        elif [ ${params.lib_type} == 3 ]; then

            printf "    third lib type\n"

        fi
        """
}


process ALIGNMENT_QC {

    tag "On: $meta.id"

    publishDir params.bismark_alignments_dir, mode: 'link'

    input:
        tuple val(meta), path(bam_file)

    output:
        tuple val(meta), path("${ meta.id }*.sorted.bam"), emit: bams
        path("${ meta.id }*_qualimap"), emit: qualimaps


    script:
    
        out_bam_file_name = params.rrbs ? "${ meta.id }_trimmed_bismark_bt2.sorted.bam" : "${ meta.id }_trimmed_bismark_bt2.deduplicated.sorted.bam"
        out_qualimap_dir = out_bam_file_name.replace(".bam", "_qualimap")

        """
        samtools sort -@ ${params.general_threads} -o ${out_bam_file_name} ${bam_file}

        qualimap bamqc -bam ${out_bam_file_name} -outdir ${out_qualimap_dir} --collect-overlap-pairs --java-mem-size=${params.qualimap_java_mem_size} -nt ${params.general_threads}
        """

}
