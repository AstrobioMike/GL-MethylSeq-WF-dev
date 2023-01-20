/*
 * Processes related to R processing with methylkit.
 */

process DIFFERENTIAL_METHYLATION_ANALYSIS {

    debug true 

    publishDir params.methylkit_outputs_dir, mode: 'link'

    input:
        file( coverage_files )
        val( simple_organism_name )
        path( runsheet )
        val( ref_org_table_link )
        val( ref_annotation_tab_link )
        val( primary_keytype )

    script:

        """

        printf "\n\n    Would be running this command:\n\n"
        
        printf "\n    Rscript --vanilla bin/differential-methylation.R \
                                        --path_to_runsheet ${ runsheet } \
                                        --simple_org_name ${ simple_organism_name } \
                                        --ref_org_table_link ${ ref_org_table_link } \
                                        --ref_annotations_tab_link ${ ref_annotation_tab_link } \
                                        --primary_keytype ${ primary_keytype }
        "

        """

}