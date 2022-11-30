/*
 * Processes for general utilities.
 */


process DOWNLOAD_GUNZIP_REFERENCES {

    // Download and decompress genome and annotation files
    tag "Organism: ${ organism_sci }  Ensembl Version: ${ensemblVersion}"
    label 'networkBound'

    publishDir params.ref_genome_dir, mode: 'link'

    input:
        tuple val(organism_sci), val(fasta_url), val(gtf_url)
        tuple val(ensemblVersion), val(ref_source)
    
    output:
        path("*.fa"), emit: fasta
        path("*.gtf"), emit: gtf

    script:
       
        // subsetting if specified
        if ( params.genomeSubsample ) {

            // doing a smaller ref genome fasta download if specific link provided
            if ( params.genomeSubsample_download_link ) {

                """
                curl -LO ${ params.genomeSubsample_download_link }
                curl -LO ${ gtf_url }

                gunzip -f *.gz

                # getting original file output names
                out_fasta=\$(ls *.fa)
                out_gtf=\$(ls *.gtf)

                # subsampling each
                samtools faidx \${out_fasta} ${ params.genomeSubsample } -o tmp.fa && mv tmp.fa \${out_fasta}

                grep '^#!' \${out_gtf} > tmp.gtf
                grep '^${ params.genomeSubsample }\t' \${out_gtf} >> tmp.gtf
                mv tmp.gtf \${out_gtf}

                # removing index
                rm *.fai
                """

            } else { 
                
                // this is the case if no specific subset fasta download link was provided, getting full one and subsetting
                """
                curl -LO ${ fasta_url }
                curl -LO ${ gtf_url }

                gunzip -f *.gz

                # getting original file output names
                out_fasta=\$(ls *.fa)
                out_gtf=\$(ls *.gtf)

                # subsampling each
                samtools faidx \${out_fasta} ${ params.genomeSubsample } -o tmp.fa && mv tmp.fa \${out_fasta}

                grep '^#!' \${out_gtf} > tmp.gtf
                grep '^${ params.genomeSubsample }\t' \${out_gtf} >> tmp.gtf
                mv tmp.gtf \${out_gtf}

                # removing index
                rm *.fai
                """

            }

        } else {

            """
            curl -LO ${ fasta_url }
            curl -LO ${ gtf_url }

            gunzip -f *.gz
            """

        }

}


process GTF_TO_PRED {

    input:
        path(ref_gtf)

    output:
        path(ref_pred), emit: pred

    script:

        // making pred outfile name
        ref_pred = ref_gtf.toString().replaceAll(".gtf", ".pred")

        """
        gtfToGenePred ${ref_gtf} ${ref_pred}
        """

}

process PRED_TO_BED {

    publishDir params.ref_genome_dir, mode: 'link'

    input:
        path(ref_pred)

    output:
        path(ref_bed), emit: bed

    script:

        // making bed outfile name
        ref_bed = ref_pred.toString().replaceAll(".pred", ".bed")

        """
        genePredToBed ${ref_pred} ${ref_bed}
        """

}

process MAKE_GENE_TRANSCRIPT_MAP {

    publishDir params.ref_genome_dir, mode: 'link'

    input:
        path(ref_gtf)

    output:
        path("*.tsv")

    script:

        // making output filename
        outfile = ref_gtf.toString().replaceAll(".gtf", "-gene-to-transcript-map.tsv")

        """
        awk ' \$3 == "transcript" ' ${ref_gtf} | cut -f 9 | tr -s ";" "\t" | \
            cut -f 1,3 | tr -s " " "\t" | cut -f 2,4 | tr -d '"' \
            > ${outfile}
        """

}
