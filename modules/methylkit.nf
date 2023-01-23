/*
 * Processes related to R processing with methylkit.
 */

process DIFFERENTIAL_METHYLATION_ANALYSIS {

    publishDir params.methylkit_outputs_dir, mode: 'link'

    input:
        path( methylkit_script )
        file( coverage_files_trigger )
        path( bismark_covs_dir )
        path( ref_dir )
        val( simple_organism_name )
        path( runsheet )
        val( ref_org_table_link )
        val( ref_annotation_tab_link )
        val( primary_keytype )

    output:
        path( "*" )

    script:

        """
        Rscript --vanilla ${ methylkit_script } \
                          --bismark_methylation_calls_dir ${ bismark_covs_dir } \
                          --path_to_runsheet ${ runsheet } \
                          --simple_org_name ${ simple_organism_name } \
                          --ref_dir ${ ref_dir } \
                          --methylkit_output_dir "." \
                          --ref_org_table_link ${ ref_org_table_link } \
                          --ref_annotations_tab_link ${ ref_annotation_tab_link } \
                          --primary_keytype ${ primary_keytype }
        """

}
