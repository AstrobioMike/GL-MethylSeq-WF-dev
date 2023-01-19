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
        path("${ params.MQCLabel }_multiqc_data.zip"), emit: zipped_report

    script:

        """
        multiqc --interactive -o ${ params.MQCLabel }_multiqc_data \
                -n ${ params.MQCLabel }_multiqc \
                --config ${ params.multiqc_config } mqc_in/*

        mv ${ params.MQCLabel }_multiqc_data/${ params.MQCLabel }_multiqc.html .
        
        zip -m -r '${ params.MQCLabel }_multiqc_data.zip' '${ params.MQCLabel }_multiqc_data'
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
        paired = params.single_end ? '' : '--paired'
        rrbs = params.rrbs ? '--rrbs' : ''
        

        """
        # this depends on the lib_type
        if [ ${ params.lib_type } == 1 ]; then

            # trimming
            trim_galore --cores ${ params.general_threads } --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

            # renaming to our convention
            if [ ${ params.single_end } == 'true' ]; then

                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz

            else

                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }*R1_trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }*R2_trimmed.fastq.gz

            fi
        
        elif [ ${ params.lib_type } == 2 ]; then

            # trimming
            trim_galore --cores ${ params.general_threads } --gzip $reads ${ non_directional } ${ rrbs } ${ paired }


            # renaming to our convention
            if [ ${ params.single_end } == 'true' ]; then

                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz

            else

                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }*R1_trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }*R2_trimmed.fastq.gz

            fi

        elif [ ${ params.lib_type } == 3 ]; then

            # the trimming command depends on if paired or not (specifically the adapter arguments)
            if [ ${ params.single_end } == 'true' ]; then

                # trimming
                trim_galore --cores ${ params.general_threads } -a AGATCGGAAGAGC --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

                # renaming to our convention
                mv ${ meta.id }*_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz                

            else

                # trimming
                trim_galore --cores ${ params.general_threads } -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC --gzip $reads ${ non_directional } ${ rrbs } ${ paired }

                # renaming to our convention
                mv ${ meta.id }*R1_trimmed.fq.gz ${ meta.id }*R1_trimmed.fastq.gz
                mv ${ meta.id }*R2_trimmed.fq.gz ${ meta.id }*R2_trimmed.fastq.gz

            fi

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
        path("${ meta.id }*_qualimap"), emit: qualimaps


    script:
    
        out_qualimap_dir = "${ bam_file }".replace("_sorted.bam", "_qualimap")

        """
        qualimap bamqc -bam ${ bam_file } -gff ${ gtf } -outdir ${ out_qualimap_dir } --collect-overlap-pairs --java-mem-size=${ params.qualimap_java_mem_size } -nt ${ params.general_threads }
        """

}
