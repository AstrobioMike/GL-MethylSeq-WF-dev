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

    publishDir params.multiqc_outputs_dir, mode: 'link', pattern: "*.zip"

    input:
        path("mqc_in/*")

    output:
        path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc.html"), emit: html
        path("${ params.MQCLabel }_multiqc_report/${ params.MQCLabel }_multiqc_data"), emit: data
        path("${ params.MQCLabel }_multiqc_report.zip"), emit: zipped_report
        path("${ params.MQCLabel }_multiqc_report"), emit: unzipped_report

    script:

        """
        multiqc --interactive -n ${ params.MQCLabel }_multiqc \
                --force --cl-config 'max_table_rows: 99999999' \
                --config ${ params.multiqc_config } \
                -o ${ params.MQCLabel }_multiqc_report mqc_in/*


        zip -r '${ params.MQCLabel }_multiqc_report.zip' '${ params.MQCLabel }_multiqc_report'
        """
}


process TRIMGALORE {

    tag "On: $meta.id"

    publishDir params.filtered_reads_dir, mode: 'link', pattern: "${ meta.id }*_trimmed.*.gz"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${ meta.id }*trimmed.*.gz"), emit: reads
        tuple val(meta), path("${ meta.id }*trimming_report.txt"), emit: reports

    script:

        non_directional = params.non_directional ? '--non_directional' : ''
        paired = "${ meta.paired_end }" == 'false' ? '' : '--paired' 

        rrbs = params.rrbs ? '--rrbs' : ''
        

        """
        # this depends on the lib_type
        if [ ${ params.lib_type } == 1 ]; then

            # trimming
            trim_galore --cores ${ params.general_threads } --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

            # renaming to our convention
            if [ ${ meta.paired_end } == 'false' ]; then

                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz

            else

                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }_R1_trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }_R2_trimmed.fastq.gz

            fi
        
        elif [ ${ params.lib_type } == 2 ]; then

            # trimming
            trim_galore --cores ${ params.general_threads } --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

            # renaming to our convention
            if [ ${ meta.paired_end } == 'false' ]; then

                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz

            else

                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }_R1_trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }_R2_trimmed.fastq.gz

            fi

        elif [ ${ params.lib_type } == 3 ]; then

            # the trimming command depends on if paired or not (specifically the adapter arguments)
            if [ ${ meta.paired_end } == 'false' ]; then

                # trimming
                trim_galore --cores ${ params.general_threads } -a AGATCGGAAGAGC --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

                # renaming to our convention
                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_pre-trimmed.fastq.gz                

            else

                # trimming
                trim_galore --cores ${ params.general_threads } -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

                # renaming to our convention
                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }_R1_pre-trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }_R2_pre-trimmed.fastq.gz

            fi

        fi
        """
}


process NUGEN_TRIM {

    tag "On: $meta.id"

    publishDir params.filtered_reads_dir, mode: 'link', pattern: "${ meta.id }*_trimmed.*.gz"

    input:
        tuple path(nugen_script), val(meta), path(reads)


    output:
        tuple val(meta), path("${ meta.id }_trimmed.*.gz"), emit: reads
        tuple val(meta), path("${ meta.id }*trim-log.txt"), emit: logs

    script:

        fastq_files = "${ meta.paired_end }" == 'false' ? "-1 ${ reads }" : "-1 ${ reads[0] } -2 ${ reads[1] }"

        """
        python ${ nugen_script } ${ fastq_files } > ${ meta.id }-NuGEN-trim-log.txt 2>&1

        # renaming output to our convention
        if [ ${ meta.paired_end } == 'false' ]; then

            mv ${ meta.id }_pre-trimmed*trimmed*.gz ${ meta.id }_trimmed.fastq.gz

        else

            mv ${ meta.id }*R1_pre-trimmed*trimmed*.gz ${ meta.id }_R1_trimmed.fastq.gz
            mv ${ meta.id }*R2_pre-trimmed*trimmed*.gz ${ meta.id }_R2_trimmed.fastq.gz

        fi
        """

}


process ALIGNMENT_QC {

    tag "On: $meta.id"

    publishDir params.bismark_alignments_dir, mode: 'link'

    input:
        tuple val(meta), path(bam_file)
        path(gtf)

    output:
        path("${ meta.id }*_qualimap.zip"), emit: qualimaps


    script:
    
        out_qualimap_dir = "${ bam_file }".replace("_sorted.bam", "_qualimap")

        """
        qualimap bamqc -bam ${ bam_file } -gff ${ gtf } -outdir ${ out_qualimap_dir } --collect-overlap-pairs --java-mem-size=${ params.qualimap_java_mem_size } -nt ${ params.general_threads }

        zip -m -r ${ out_qualimap_dir }.zip ${ out_qualimap_dir }
        """

}
