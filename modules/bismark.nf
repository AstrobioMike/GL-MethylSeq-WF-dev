/*
 * Processes related to bismark processing.
 */

process GEN_BISMARK_REF {

    publishDir "./", mode: 'link', pattern: "*.zip"

    input:
        path( input_ref_fasta )

    output:
        path( params.bismark_index_dir.toString().replaceAll("/", ".zip") )
        path( params.bismark_index_dir ), emit: ch_bismark_index_dir

    script:

        zip_out = params.bismark_index_dir.toString().replaceAll("/", ".zip")

        """
        mkdir -p ${ params.bismark_index_dir }

        # copying ref fasta in there so everything exists with the bismark index dir (required by bismark)
        cp ${ input_ref_fasta } ${ params.bismark_index_dir }

        # making index
        bismark_genome_preparation --bowtie2 --parallel ${ params.bismark_index_creation_threads } ${ params.bismark_index_dir } > bismark-genome-preparation-output.txt 2>&1

        # moving log file into bismark dir since i can't figure out how to do it with publishDirs above
        mv bismark-genome-preparation-output.txt ${ params.bismark_index_dir }

        # creating genomic (di)-nucleotide stats file for later use by alignment step with --nucleotide_coverage flags
        bam2nuc --genome_folder ${ params.bismark_index_dir } --genomic_composition_only

        # zipping up for packaging with all outputs
        zip -r ${ zip_out } ${ params.bismark_index_dir }
        """        

}


process ALIGN {

    tag "On: $meta.id"

    publishDir params.bismark_alignments_dir, mode: 'link', pattern: "${ meta.id }*_sorted.bam"
    publishDir params.bismark_alignments_dir, mode: 'link', pattern: "${ meta.id }*.nucleotide_stats.txt"


    input:
        tuple val(meta), path(reads), path( bismark_index_dir )

    output:
        tuple val(meta), path("${ meta.id }_bismark_*_sorted.bam"), emit: bams
        tuple val(meta), path("${ meta.id }_bismark_*_report.txt"), emit: reports
        tuple val(meta), path("${ meta.id }_bismark_*.nucleotide_stats.txt"), emit: nuc_stats

    script:

        non_directional = params.non_directional ? '--non_directional' : ''
        fastq_files = "${ meta.paired_end }" == 'false' ? "${ reads }" : "-1 ${ reads[0] } -2 ${ reads[1] }" 


        """
        bismark --bam --non_bs_mm --nucleotide_coverage --parallel ${ params.bismark_align_threads } -p ${ params.bismark_align_threads } --genome_folder ${ bismark_index_dir } ${ non_directional } ${ fastq_files } 

        sorted_bam_name=\$(ls *.bam | sed 's/_trimmed//' | sed 's/.bam/_sorted.bam/')

        # sorting bam file
        samtools sort -@ ${ params.general_threads } -o \${sorted_bam_name} *.bam

        # renaming report files so it is cleaner and match up with autodetection of sorted bam files by bismark2summary later
        new_report_filename=\$(ls *_report.txt | sed 's/_trimmed//' | sed 's/_\\(.E\\)/_sorted_\\1/')
        mv *_report.txt \${new_report_filename}

        # renaming nucleotide stats file to be cleaner
        new_stats_filename=\$(ls *.nucleotide_stats.txt | sed 's/_trimmed//')
        mv *.nucleotide_stats.txt \${new_stats_filename}
        """

}


process DEDUPLICATE {

    tag "On: $meta.id"

    publishDir params.bismark_alignments_dir, mode: 'link', pattern: "${ meta.id }_bismark_*.bam"

    input:
        tuple val(meta), path(bam_file)

    output:
        tuple val(meta), path("${ meta.id }_bismark_*.bam"), emit: bams
        tuple val(meta), path("${ meta.id }*_report.txt"), emit: reports

    script:

        """
        deduplicate_bismark ${ bam_file }
        """

}


process EXTRACT_METHYLATION_CALLS {

    tag "On: $meta.id"

    publishDir params.bismark_methylation_calls_dir, mode: 'link', pattern: "*cov.gz"
    publishDir params.bismark_methylation_calls_dir, mode: 'link', pattern: "*.zip"

    input:
        tuple val(meta), path(bam_file)
        path(bismark_index_dir)


    output:
        tuple val(meta), path("${ meta.id }*.cov.gz"), emit: covs
        tuple val(meta), path("${ meta.id }*.bedGraph.gz"), emit: beds
        tuple val(meta), path("*${ meta.id }*.txt.gz"), emit: contexts
        tuple val(meta), path("${ meta.id }*.M-bias.txt"), emit: biases
        tuple val(meta), path("${ meta.id }*_report.txt"), emit: reports
        tuple val(meta), path("${ meta.id }*cytosine_context_summary.txt"), emit: context_summary
        tuple val(meta), path("${ meta.id }*.zip"), emit: zip

    script:

        additional_args = "${ meta.paired_end }" == 'false' ? "" : "--ignore_r2 2 --ignore_3prime_r2 2"
        zip_out = "${ meta.id }_bismark_methylation_calls.zip"


        """
        bismark_methylation_extractor --bedGraph --gzip --comprehensive --cytosine_report --genome_folder ${ bismark_index_dir } ${ additional_args } ${ bam_file }
        zip ${ zip_out } *.gz *.txt
        """

}


process GEN_BISMARK_SAMPLE_REPORT {

    tag "On: $meta.id"

    publishDir params.individual_sample_reports, mode: 'link', pattern: "*.html"

    input:
        tuple val(meta), path(alignment_report), path(meth_calls_report), path(m_bias_report), path(nucleotide_covs), file(dedupe_report)

    output:
        tuple val(meta), path("${ meta.id }*_report.html"), emit: reports

    script:

        additional_args = params.rrbs ? "" : "--dedup_report ${ dedupe_report }"

        """
        bismark2report --alignment_report ${ alignment_report } --splitting_report ${ meth_calls_report } --mbias_report ${ m_bias_report } --nucleotide_report ${ nucleotide_covs } ${ additional_args }
        """

}


process GEN_BISMARK_SUMMARY {

    publishDir params.bismark_summary_dir, mode: 'link', pattern: "bismark_summary_report.*"

    input:
        file(all_bams_and_reports)

    output:
        tuple path("bismark_summary_report.html"), path("bismark_summary_report.txt"), emit: reports

    script:

        """
        bismark2summary *.bam
        """

}
