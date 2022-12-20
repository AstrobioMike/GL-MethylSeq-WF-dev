library(optparse)
## other packages loaded below after help menu might be called

############################################
############ handling arguments ############
############################################

parser <- OptionParser(description = "\n  This is really only intented to be called from within the GeneLab MethylSeq workflow.")

parser <- add_option(parser, c("-v", "--verbose"),
                     action = "store_true",
                     default = FALSE, help = "Print extra output")

parser <- add_option(parser, c("--bismark_methylation_calls_dir"),
                     default = "Bismark_Methylation_calls",
                     help = "Directory holding *bismark.cov.gz files")

parser <- add_option(parser, c("--metadata_dir"),
                     default = "Metadata",
                     help = "Directory holding metadata files")

parser <- add_option(parser, c("--ref_dir"),
                     default = "Reference_Genome_Files",
                     help = "Directory holding reference genome files (e.g. *.bed and *.gtf files)")

parser <- add_option(parser, c("--methylkit_output_dir"),
                     default = "MethylKit_Outputs",
                     help = "Directory for methylkit output files")

parser <- add_option(parser, c("--limit_samples_to"),
                     default = "all",
                     help = "Limits the number of samples being processed (won't do real factor comparisons if set)")

parser <- add_option(parser, c("--ref_genome_string"),
                     help = "Reference genome used (just for recording, not used here)")

parser <- add_option(parser, c("--ref_annotations_tab_link"),
                     help = "Link to reference-genome annotations table")

parser <- add_option(parser, c("--mc_cores"), default = 4, type = "integer",
                     help = "Passed to mc.cores argument of calculateDiffMeth() call")

parser <- add_option(parser, c("--getMethylDiff_difference"), default = 25, type = "integer",
                     help = "Passed to difference argument of getMethylDiff() call")

parser <- add_option(parser, c("--getMethylDiff_qvalue"), default = 0.01, type = "double",
                     help = "Passed to qvalue argument of getMethylDiff() call")

parser <- add_option(parser, c("--primary_keytype"),
                     help = "The keytype to use for mapping annotations (usually 'ENSEMBL' for most things; 'TAIR' for plants)")

args <- parse_args(parser)

############################################
########### for testing purposes ###########
############################################

# args$v <- TRUE
# args$bismark_methylation_calls_dir <- "Bismark_Methylation_Calls"
# args$metadata_dir <- "Metadata"
# args$ref_dir <- "Reference_Genome_Files"
# args$methylkit_output_dir <- "MethylKit_Outputs"
# args$limit_samples_to <- 6
# args$ref_genome_string <- "Mmus_GRCm39"
# args$ref_annotations_tab_link <- "https://figshare.com/ndownloader/files/36597114"
# args$getMethylDiff_difference <- 1
# args$getMethylDiff_qvalue <- 0.9
# args$primary_keytype <- "ENSEMBL"

############################################


############################################
########### checking on arguments ##########
############################################

# checking required arguments were set
required_args <- c("ref_genome_string" = "--ref_genome_string",
                   "ref_annotations_tab_link" = "--ref_annotations_tab_link",
                   "primary_keytype" = "--primary_keytype")


for ( arg in names(required_args) ) {

    tryCatch( { get(arg, args) }, error = function(e) { 
              
                      stop(cat("\nThe '", as.character(required_args[arg]),
                               "' argument must be provided. Cannot proceed.\n", sep = ""), call. = FALSE)

    })
    
}

# checking primary keytype is what's expected
currently_accepted_keytypes <- c("ENSEMBL", "TAIR")
if ( ! args$primary_keytype %in% currently_accepted_keytypes ) {
    
    stop(cat("\nThe current potential --primary_keytypes are:", currently_accepted_keytypes,
             "\nCannot proceed with:", args$primary_keytype, "\n"), call. = FALSE)
}

############################################
############# loading packages #############
############################################

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(methylKit)))
suppressWarnings(suppressMessages(library(genomation)))

############################################


############################################
############# helper functions #############
############################################
get_single_file_path <- function(target_dir, search_pattern) {
    
    hits <- list.files(target_dir, pattern = search_pattern, full.names = TRUE)
    
    # checking only one was found matching pattern search
    if ( length(hits) != 1 ) { 
        
        error_message <- cat("\nA single file was not found in the ", target_dir,
                             " directory based on the search pattern '",
                             search_pattern, "'. Cannot proceed.\n", sep = "")
        
        stop(error_message, call. = FALSE)
    }
    
    return(hits)
    
}

order_input_files <- function(sample_names, paths) {
    
    # this function takes in a vector of sample_names and
    # a vector of paths to coverage files, and returns a vector
    # of coverage files in the same order as the sample names
    
    ordered_paths <- c()
    
    for ( sample in sample_names ) {
        
        search_pattern <- paste0(sample, "_trimmed")
        hits <- paths[grep(search_pattern, paths)]
        
        # making sure there is exactly one match
        if ( length(hits) != 1 ) {
            
            stop("\nThere was a problem matching up sample names with their coverage files. Cannot proceed.\n", call. = FALSE)
            
        }
        
        ordered_paths <- c(ordered_paths, hits)
        
    }
    
    return(ordered_paths)
    
}

############################################


############################################
############# setting things up ############
############################################

### finding reference bed file
ref_bed_path <- get_single_file_path(args$ref_dir, ".*.bed")


### finding reference gene-to-transcript mapping file
ref_gene_transcript_map_path <- get_single_file_path(args$ref_dir, ".*-gene-to-transcript-map.tsv")


### getting path to runsheet
runsheet_path <- get_single_file_path(args$metadata_dir, ".*_runsheet.csv")

# reading runsheet
runsheet <- read.csv(runsheet_path)

### getting all factors
# mock factors if wanted for testing
# runsheet$Factor.Value.Other <- c(rep("Mad", 10), rep("Dog", 6))
# runsheet$Factor.Value.Other2 <- c(rep("LilMac", 6), rep("Fighter", 10))
# runsheet$Factor.Value.Other2 <- c(rep("A", 6), rep("B", 8), rep("C", 2))
factors <- runsheet %>% dplyr::select(starts_with("Factor.Value"))
colnames(factors) = paste("factor", 1:dim(factors)[2], sep = "_")

# checking if there is more than two unique values in a given factor, as methylkit isn't built for that
    # might not need this after doing combined only way as RNAseq does, need to check later
for ( factor in colnames(factors)) {

    curr_unique_entries <- factors[factor] %>% unique() %>% pull()

    if ( length(curr_unique_entries) > 2 ) {

        error_message <- cat("\nOne of the factors has more than two entries:\n\n  ",
                             print(curr_unique_entries),
                             "\n\nOnly pairwise comparisons can be done currently. Exiting for now.\n")
        stop(error_message, call. = FALSE)

    }

}


# factor dataframe
factor_df <- data.frame(sample_id = runsheet %>% pull("Sample.Name"), factors)

# subsetting runsheet if specified
if ( args$limit_samples_to != "all" ) {
    
    args$limit_samples_to <- as.integer(args$limit_samples_to)

    runsheet <- runsheet[ 1:args$limit_samples_to, ]
    
    # making up factors for testing if things were subset (so that we have multiple factor types even if not really in the subset samples)
    len_type_A <- ceiling(args$limit_samples_to / 2)
    len_type_B <- args$limit_samples_to - len_type_A
    
    factor_1 <- c(rep("A", len_type_A), rep("B", len_type_B))
    
    factor_df <- data.frame(sample_id = runsheet %>% pull("Sample.Name"), factor_1)

}

if ( args$v ) { cat("\nFactor df:\n") ; print(factor_df) ; cat("\n") }

# storing just sample names in a vector
sample_names <- runsheet$Sample.Name

# getting list of *.bismark.cov.gz files
bismark_cov_paths <- list.files(args$bismark_methylation_calls_dir, pattern = ".*.bismark.cov.gz", full.names = TRUE)

# making sure files list matches length of runsheet
if ( dim(runsheet)[1] != length(bismark_cov_paths) ) {

    error_message <- cat("\nThe number of '*.bismark.cov.gz' files found in the", 
                         args$bismark_methylation_calls_dir, 
                         "directory does not match the number of samples specified in the runsheet. Cannot proceed.\n")
    stop(error_message, call. = FALSE)
    
}

if ( args$v ) { cat("\nBismark coverage files:\n") ; print(bismark_cov_paths) ; cat("\n") }

# making sure coverage-file-paths vector is in the same order as the sample names
bismark_cov_paths <- order_input_files(sample_names, bismark_cov_paths)

# making table with filenames and bismark coverage file path
samples_and_covs_paths_df <- data.frame("sample_id" = sample_names, "coverage_file_path" = bismark_cov_paths)

## GeneLab combines all factors and just runs those contrasts, so doing things that way
# a lot of this initial structuring comes from Jonathan Oribello's work in the GeneLab RNAseq workflow
study_df <- factor_df %>% column_to_rownames("sample_id")

# if there are multiple factors, here we are concatenating them to make one combined one
groups <- apply(study_df, 1, paste, collapse = " & ")
group_names <- paste0("(", groups, ")", sep = "") # human readable group names
safe_group_names <- make.names(groups) # group naming compatible with R models
names(safe_group_names) <- group_names

samples_and_combined_factors_df <- data.frame("sample_id" = sample_names, "combined_factor" = safe_group_names)

if ( args$v ) { cat("\nSamples and factors df:\n") ; print(samples_and_combined_factors_df) ; cat("\n") }


##### Format contrasts table, defining pairwise comparisons for all groups #####
contrasts <- combn(levels(factor(safe_group_names)), 2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(safe_group_names))), 2)


####### NOTE TO MIKE DURING DEV #######
# this is about output A vs B AND B vs A, as is currently done with RNAseq (both are done and output)
    # unlike with deseq2, each contrast needs to be re-run to have the outputs the other way
    # not all the output info can go in one table (by just adding columns), so we are doubling all files by reporting both
# so i think ultimately we shouldn't do all contrasts both ways with methylseq data
#######################################

## this way would be if doing A vs B and B vs A
# contrast.names <- c(paste(contrast.names[1,], contrast.names[2,], sep = "v"), 
#                     paste(contrast.names[2,], contrast.names[1,], sep = "v"))
# contrasts <- cbind(contrasts, contrasts[c(2,1),])

## this way is only doing one-way contrasts
contrast.names <- paste(contrast.names[1,], contrast.names[2,], sep = "v")
colnames(contrasts) <- contrast.names

if ( args$v ) { cat("\nContrasts table:\n") ; print(contrasts) ; cat("\n") }

# making a single table with info needed for methylkit
sample_meth_info_df <- full_join(samples_and_covs_paths_df, samples_and_combined_factors_df, "sample_id")

if ( args$v ) { cat("\nPrimary methylkit df:\n") ; print(sample_meth_info_df) ; cat("\n") }

# making output directory
dir.create(args$methylkit_output_dir, showWarnings = FALSE)

# reading in transcript features
gene.obj <- readTranscriptFeatures(ref_bed_path, up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, unique.prom = TRUE)

# reading in gene to transcript mapping file
gene_transcript_map <- 
    read.table(ref_gene_transcript_map_path, sep = "\t", col.names = c("gene_ID", "feature.name"))

# reading in functional annotation table
options(timeout = 600)

functional_annots_tab <- 
    read.table(args$ref_annotations_tab_link, sep = "\t", quote = "", header = TRUE)

############################################


############################################
########## moving onto methylkit ##########
############################################


### looping through contrasts and running methylkit and creating outputs for each ###
for ( i in 1:dim(contrasts)[2]) { 

    # getting current contrast
    curr_contrasts_vec <- contrasts[, i]
    
    # prefix for output files
    curr_output_prefix <- paste(curr_contrasts_vec, collapse = "_vs_")
    
    # making directory for output files of current contrast
    curr_output_dir <- file.path(args$methylkit_output_dir, curr_output_prefix)
    dir.create(curr_output_dir, showWarnings = FALSE)
    
    # getting subset sample info table relevant to current contrast
    curr_sample_info_df <- sample_meth_info_df %>% filter(combined_factor %in% curr_contrasts_vec)

    # getting which samples are relevant to current contrast
    curr_samples_vec <- curr_sample_info_df %>% pull(sample_id)
    
    # getting which files are relevant to current contrast
    curr_coverage_file_paths <- curr_sample_info_df %>% pull(coverage_file_path)
    
    # making binary vector for treatment argument to methRead()
    curr_treatment_vec <- c()
    for ( value in curr_sample_info_df$combined_factor ) {
        
        if ( value == curr_sample_info_df$combined_factor[1] ) {
            
            curr_treatment_vec <- c(curr_treatment_vec, 1)
            
        } else {
            
            curr_treatment_vec <- c(curr_treatment_vec, 0)
            
        }

    }
    

    ### setting up methylkit object
    curr_obj <- methRead(location = as.list(curr_coverage_file_paths),
                         sample.id = as.list(curr_samples_vec),
                         treatment = curr_treatment_vec,
                         pipeline = "bismarkCoverage",
                         assembly = args$ref_genome_string,
                         header = FALSE,
                         mincov = 10)
    
    ### Individual-base analysis
    # merging samples
    curr_meth <- unite(curr_obj)
    
    # calculating differential methylation
    curr_myDiff <- calculateDiffMeth(curr_meth, mc.cores = args$mc_cores)
    
    # getting just sig hyper-methlated bases
    curr_myDiff25p.hyper <- getMethylDiff(curr_myDiff, difference = args$getMethylDiff_difference, 
                                          qvalue = args$getMethylDiff_qvalue, type = "hyper")

    # getting df of just sig hyper-methylated bases
    curr_sig_bases_hyper_tab <- getData(curr_myDiff25p.hyper) %>% arrange(qvalue)
    
    # getting just sig hypo-methylated bases
    curr_myDiff25p.hypo <- getMethylDiff(curr_myDiff, difference = args$getMethylDiff_difference,
                                         qvalue = args$getMethylDiff_qvalue, type = "hypo")

    # getting df of just sig hypo-methylated bases
    curr_sig_bases_hypo_tab <- getData(curr_myDiff25p.hypo) %>% arrange(qvalue)
    
    # getting all sig differentially methylated bases
    curr_myDiff25p <- getMethylDiff(curr_myDiff, difference = args$getMethylDiff_difference,
                                    qvalue = args$getMethylDiff_qvalue)

    # getting df of all sig methylated bases
    curr_sig_bases_all_tab <- getData(curr_myDiff25p) %>% arrange(qvalue)

    
    ### Tile analysis ###
    # tiling
    curr_tiles_obj <- tileMethylCounts(curr_obj, win.size = 1000, step.size = 1000, cov.bases = 10)
    
    # merging tiled samples
    curr_tiles_meth <- unite(curr_tiles_obj)
    
    # calculating differential methylation on tiles
    curr_tiles_diff <- calculateDiffMeth(curr_tiles_meth, mc.cores = args$mc_cores)
    
    # getting sig hyper-methylated tiles
    curr_tiles_myDiff25p.hyper <- getMethylDiff(curr_tiles_diff, difference = args$getMethylDiff_difference, 
                                                qvalue = args$getMethylDiff_qvalue, type = "hyper")

    # getting table of sig hyper-methylated tiles
    curr_tiles_sig_hyper_tab <- getData(curr_tiles_myDiff25p.hyper) %>% arrange(qvalue)
    
    # getting sig hypo-methylated tiles
    curr_tiles_myDiff25p.hypo <- getMethylDiff(curr_tiles_diff, difference = args$getMethylDiff_difference, 
                                               qvalue = args$getMethylDiff_qvalue, type = "hypo")

    # getting table of sig hypo-methylated tiles
    curr_tiles_sig_hypo_tab <- getData(curr_tiles_myDiff25p.hypo) %>% arrange(qvalue)
    
    # getting all sig differentially methylated tiles
    curr_tiles_myDiff25p <- getMethylDiff(curr_tiles_diff, difference = args$getMethylDiff_difference, 
                                          qvalue = args$getMethylDiff_qvalue)

    # making table of all sig differentially methylated tiles
    curr_sig_tiles_all_tab <- getData(curr_tiles_myDiff25p) %>% arrange(qvalue)
    
    
    ### Adding feature information ###
    
    ## adding features to individual-base objects
    curr_diffAnn <- annotateWithGeneParts(as(curr_myDiff25p, "GRanges"), gene.obj)
    curr_diffAnn.hyper <- annotateWithGeneParts(as(curr_myDiff25p.hyper, "GRanges"), gene.obj)
    curr_diffAnn.hypo <- annotateWithGeneParts(as(curr_myDiff25p.hypo, "GRanges"), gene.obj)
    
    # making base-level sig tables with features 
    curr_sig_all_bases_tab_with_features <- cbind(data.frame(curr_myDiff25p), 
                                                  getAssociationWithTSS(curr_diffAnn), 
                                                  as.data.frame(getMembers(curr_diffAnn))) %>% .[,-c(8)]

    curr_sig_hyper_bases_tab_with_features <- cbind(data.frame(curr_myDiff25p.hyper), 
                                                    getAssociationWithTSS(curr_diffAnn.hyper), 
                                                    as.data.frame(getMembers(curr_diffAnn.hyper))) %>% .[,-c(8)]
    
    curr_sig_hypo_bases_tab_with_features <- cbind(data.frame(curr_myDiff25p.hypo), 
                                                   getAssociationWithTSS(curr_diffAnn.hypo), 
                                                   as.data.frame(getMembers(curr_diffAnn.hypo))) %>% .[,-c(8)]
    
    
    ## adding features to tiles objects
    curr_tiles_diffAnn <- annotateWithGeneParts(as(curr_tiles_myDiff25p, "GRanges"), gene.obj)
    curr_tiles_diffAnn.hyper <- annotateWithGeneParts(as(curr_tiles_myDiff25p.hyper, "GRanges"), gene.obj)
    curr_tiles_diffAnn.hypo <- annotateWithGeneParts(as(curr_tiles_myDiff25p.hypo, "GRanges"), gene.obj)
    
    # making tiles sig tables with features 
    curr_tiles_sig_all_tab_with_features <- cbind(data.frame(curr_tiles_myDiff25p), 
                                                  getAssociationWithTSS(curr_tiles_diffAnn), 
                                                  as.data.frame(getMembers(curr_tiles_diffAnn))) %>% .[,-c(8)]
    
    curr_tiles_sig_hyper_tab_with_features <- cbind(data.frame(curr_tiles_myDiff25p.hyper), 
                                                    getAssociationWithTSS(curr_tiles_diffAnn.hyper), 
                                                    as.data.frame(getMembers(curr_tiles_diffAnn.hyper))) %>% .[,-c(8)]
    
    curr_tiles_sig_hypo_tab_with_features <- cbind(data.frame(curr_tiles_myDiff25p.hypo), 
                                                   getAssociationWithTSS(curr_tiles_diffAnn.hypo), 
                                                   as.data.frame(getMembers(curr_tiles_diffAnn.hypo))) %>% .[,-c(8)]

    
    ### Adding functional annotations ###
    
    ## for individual-base outputs
    # for each transcript ID in the sig_all_bases_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table
    curr_sig_all_bases_tab_with_features_and_gene_IDs <- 
        left_join(curr_sig_all_bases_tab_with_features, gene_transcript_map)

    curr_sig_hyper_bases_tab_with_features_and_gene_IDs <- 
        left_join(curr_sig_hyper_bases_tab_with_features, gene_transcript_map)

    curr_sig_hypo_bases_tab_with_features_and_gene_IDs <- 
        left_join(curr_sig_hypo_bases_tab_with_features, gene_transcript_map)
    
    # now adding full annotations
    curr_sig_all_bases_tab_with_features_and_annots <- 
        left_join(curr_sig_all_bases_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = args$primary_keytype))

    curr_sig_hyper_bases_tab_with_features_and_annots <- 
        left_join(curr_sig_hyper_bases_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = args$primary_keytype))
    
    curr_sig_hypo_bases_tab_with_features_and_annots <- 
        left_join(curr_sig_hypo_bases_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = args$primary_keytype))
    
            
    # and writing out
    curr_sig_all_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-bases.tsv"))
    
    write.table(curr_sig_all_bases_tab_with_features_and_annots, curr_sig_all_bases_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    curr_sig_hyper_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-bases.tsv"))
    
    write.table(curr_sig_hyper_bases_tab_with_features_and_annots, curr_sig_hyper_bases_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    curr_sig_hypo_bases_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-bases.tsv"))
    
    write.table(curr_sig_hypo_bases_tab_with_features_and_annots, curr_sig_hypo_bases_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    ## for tiles output
    # for each transcript ID in the tiles_sig_all_out_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table
    curr_sig_all_tiles_tab_with_features_and_gene_IDs <- 
        left_join(curr_tiles_sig_all_tab_with_features, gene_transcript_map)
    
    curr_sig_hyper_tiles_tab_with_features_and_gene_IDs <- 
        left_join(curr_tiles_sig_hyper_tab_with_features, gene_transcript_map)
    
    curr_sig_hypo_tiles_tab_with_features_and_gene_IDs <- 
        left_join(curr_tiles_sig_hypo_tab_with_features, gene_transcript_map)
    
    # now adding full annotations
    curr_sig_all_tiles_tab_with_features_and_annots <- 
        left_join(curr_sig_all_tiles_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
    
    curr_sig_hyper_tiles_tab_with_features_and_annots <- 
        left_join(curr_sig_hyper_tiles_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
    
    curr_sig_hypo_tiles_tab_with_features_and_annots <- 
        left_join(curr_sig_hypo_tiles_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
    
    # and writing out
    curr_sig_all_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-tiles.tsv"))
    
    write.table(curr_sig_all_tiles_tab_with_features_and_annots, curr_sig_all_tiles_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    curr_sig_hyper_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-tiles.tsv"))
    
    write.table(curr_sig_hyper_tiles_tab_with_features_and_annots, curr_sig_hyper_tiles_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    curr_sig_hypo_tiles_tab_with_features_and_annots_path <- 
        file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-tiles.tsv"))
    
    write.table(curr_sig_hypo_tiles_tab_with_features_and_annots, curr_sig_hypo_tiles_tab_with_features_and_annots_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    ### Making overview figure of percent diff. methylation across features ###
    curr_sig_diff_bases_across_features_plot_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-bases-across-features.pdf"))
    pdf(curr_sig_diff_bases_across_features_plot_path)
    plotTargetAnnotation(curr_diffAnn, precedence = TRUE, main = "% of sig. diff. methylated sites across features")
    dev.off()
    
    curr_sig_diff_tiles_across_features_plot_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-tiles-across-features.pdf"))
    pdf(curr_sig_diff_tiles_across_features_plot_path)
    plotTargetAnnotation(curr_tiles_diffAnn, precedence = TRUE, main = "% of sig. diff. methylated tiles across features")
    dev.off()
    
}

### making and writing out a table of base-level percent methylated (treatment vector doesn't matter here, just making a mock one)
len_1s <- ceiling(length(sample_names) / 2)
len_0s <- length(sample_names) - len_1s
mock_treatment_vec <- c(rep(1, len_1s), rep(0, len_0s))

obj <- methRead(location = as.list(sample_meth_info_df %>% pull(coverage_file_path)),
                sample.id = as.list(sample_meth_info_df %>% pull(sample_id)),
                treatment = mock_treatment_vec,
                pipeline = "bismarkCoverage",
                assembly = args$ref_genome_string,
                header = FALSE,
                mincov = 10)

meth <- unite(obj)

perc.meth <- percMethylation(meth, rowids = TRUE)
perc.meth <- perc.meth %>% data.frame() %>% rownames_to_column("location")

# writing out
perc.meth_path <- file.path(args$methylkit_output_dir, "base-level-percent-methylated.tsv")
write.table(perc.meth, perc.meth_path, sep = "\t", quote = FALSE, row.names = FALSE)

### making and writing out a table of tile-level percent methylated (contrasts don't matter here)
tiles_obj <- tileMethylCounts(obj, win.size = 1000, step.size = 1000, cov.bases = 10)
tiles_meth <- unite(tiles_obj)

tiles_perc.meth <- percMethylation(tiles_meth, rowids = TRUE)
tiles_perc.meth <- tiles_perc.meth %>% data.frame() %>% rownames_to_column("location")

tiles_perc.meth_path <- file.path(args$methylkit_output_dir, "tile-level-percent-methylated.tsv")
write.table(tiles_perc.meth, tiles_perc.meth_path, sep = "\t", quote = FALSE, row.names = FALSE)
