/* Processes dealing with retrieving data from GeneLab */
/* Modified from Jonathan Oribello: https://github.com/asaravia-butler/GeneLab_Data_Processing/blob/b6576eade4b77e3bbf5bcad573c8d77fa7d20d1c/RNAseq/Workflow_Documentation/NF_RCP-F/workflow_code/NF_RCP-F_1.0.0/modules/genelab.nf */

// TODO Migrate these CLI args to updated dp tools package API
process METHYLSEQ_RUNSHEET_FROM_GLDS {

    // Downloads isazip and creates run sheets using GeneLab API
    tag "${ glds_accession }"
    publishDir params.metadata_dir, mode: "link"

    input:
        val(glds_accession)

    output:
        path("${ glds_accession }*.csv"), emit: runsheet
        path("*.zip"), emit: isazip

    script:
        """
        dpt-get-isa-archive --accession ${ glds_accession }
        dpt-isa-to-runsheet --accession ${ glds_accession } \
            --config-type methylSeq --isa-archive *.zip
        """
}


process STAGE_RAW_READS {

    // Stages the raw reads into appropriate publish directory
    tag "${ meta.id }"

    input:
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        if ( meta.paired_end ) {
            """
            cp -L 1.gz ${meta.id}_R1_raw.fastq.gz
            cp -L 2.gz ${meta.id}_R2_raw.fastq.gz
            """
        } else {
            """
            cp -L 1.gz  ${meta.id}_raw.fastq.gz
            """
        }

}


process GENERATE_METASHEET {

    // Generates a metadata table, not used in further processing
    tag "${ params.gldsAccession }"
    publishDir params.metadata_dir, mode: "link"

    input:
        path("isa.zip")
        path(runsheet)

    output:
        path("${ params.gldsAccession }_metadata_table.tsv"), emit: metasheet

    script:
        """
        create_table_v2.py --accession ${ params.gldsAccession }  \
                        --isa-zip isa.zip \
                        --output-dir . \
                        --runsheet ${ runsheet }
        """

}


process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    tuple val(organism_sci), val(fasta_url), val(gtf_url), emit: reference_genome_urls
    val(annotations_db_url), emit: annotations_db_url
    tuple val(ensemblVersion), val(ensemblSource), emit: reference_version_and_source
    val(simple_organism_name), emit: simple_organism_name
  
  exec:
    def organisms = [:]
    // println "  ${YELLOW}Fetching table from ${annotations_csv_url_string}${NC}"
    
    // download data to memory
    annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    fasta_url = organisms[organism_key][5]
    gtf_url = organisms[organism_key][6]
    annotations_db_url = organisms[organism_key][9]
    ensemblVersion = organisms[organism_key][3]
    ensemblSource = organisms[organism_key][4]
    simple_organism_name = organisms[organism_key][0]

    // println "  ${YELLOW}PARSING_ANNOTATIONS_TABLE:"
    // println "    Values parsed for '${organism_key}':"
    // println "--------------------------------------------------"
    // println "      - fasta_url: ${fasta_url}"
    // println "      - gtf_url: ${gtf_url}"
    // println "      - annotations_db_url: ${annotations_db_url}"
    // println "      - ensemblVersion: ${ensemblVersion}"
    // println "      - ensemblSource: ${ensemblSource}"
    // println "--------------------------------------------------${NC}"

}

// process GENERATE_MD5SUMS {
//   // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
//   tag "${ params.gldsAccession }"
//   publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
//     mode: params.publish_dir_mode

//   input:
//     path(data_dir)
//     path(runsheet)

//   output:
//     path("*md5sum*")
//     path("Missing_md5sum_files.txt"), optional: true

//   script:
//     """
//     generate_md5sum_files.py  --root-path ${ data_dir } \\
//                               --runsheet-path ${ runsheet }
//     """
// }

// process UPDATE_ISA_TABLES {
//   // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
//   tag "${ params.gldsAccession }"
//   publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
//     mode: params.publish_dir_mode

//   input:
//     path(data_dir)
//     path(runsheet)

//   output:
//     path("updated_curation_tables") // directory containing extended ISA tables

//   script:
//     """
//     update_curation_table.py  --root-path ${ data_dir } \\
//                               --runsheet-path ${ runsheet }
//     """
// }

// process SOFTWARE_VERSIONS {
//   // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
//   tag "${ params.gldsAccession }"
//   publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
//     mode: params.publish_dir_mode

//   input:
//     path("software_versions.txt")

//   output:
//     path("software_versions.md")

//   script:
//     """
//     format_software_versions.py software_versions.txt
//     """
// }

// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Original Function Credit: Dr. Harshil Patel
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def get_runsheet_paths(LinkedHashMap row) {

    def ORGANISMS = ["mus_musculus":"MOUSE",
                     "danio_rerio":"ZEBRAFISH",
                     "rattus_norvegicus":"RAT",
                     "homo_sapiens":"HUMAN",
                     "drosophila_melanogaster":"FLY",
                     "caenorhabditis_elegans":"WORM",
                     "arabidopsis_thaliana":"ARABIDOPSIS"]

    def PRIMARY_KEYS = ["mus_musculus":"ENSEMBL",
                        "danio_rerio":"ENSEMBL",
                        "rattus_norvegicus":"ENSEMBL",
                        "homo_sapiens":"ENSEMBL",
                        "drosophila_melanogaster":"ENSEMBL",
                        "caenorhabditis_elegans":"ENSEMBL",
                        "arabidopsis_thaliana":"TAIR"]

    def meta = [:]
    meta.id                         = row["Sample Name"]
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.primary_keytype            = PRIMARY_KEYS[meta.organism_sci]
    meta.paired_end                 = row.paired_end.toBoolean()

    def array = []
    def raw_reads = []


    raw_reads.add(file(row.read1_path))

    if (meta.paired_end) {

        raw_reads.add(file(row.read2_path))

    }
    
    array = [meta, raw_reads]

    return array

}
