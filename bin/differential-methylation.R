library(optparse)
library(tidyverse)
library(methylKit)
library(genomation)

parser <- OptionParser()

parser <- add_option(parser, c("-v", "--verbose"),
                     action = "store_true",
                     default = FALSE, help = "Print extra output")


# directory of Bismark_Methylation_Calls
parser <- add_option(parser, c("--bismark_methylation_calls_dir"),
                     default = "Bismark_Methylation_calls",
                     help = "Directory holding *bismark.cov.gz files")

# directory of metadata
parser <- add_option(parser, c("--metadata_dir"),
                     default = "Metadata",
                     help = "Directory holding metadata files")


# directory of reference genome files
parser <- add_option(parser, c("--ref_dir"),
                     default = "Reference_Genome_Files",
                     help = "Directory holding reference genome files (e.g. *.bed and *.gtf files)")


# methylkit output directory
parser <- add_option(parser, c("--methylkit_output_dir"),
                     default = "MethylKit_Outputs",
                     help = "Directory holding metadata files")


# directory of metadata
parser <- add_option(parser, c("--limit_samples_to"),
                     default = "all",
                     help = "Limits the number of samples being processed (won't do real factor comparisons if set)")


# reference genome used (just for recording in methylkit)
parser <- add_option(parser, c("--ref_genome_string"),
                     help = "Reference genome used (just for recording, not used here)")


# number of threads to run in parallel for methylkit calls
parser <- add_option(parser, c("--mc_cores"), default = 4, type = "integer",
                     help = "Passed to mc.cores argument of calculateDiffMeth() call")


args <- parse_args(parser)

# checking required arguments were set
if ( is.na(args$ref_genome_string )) {
    stop("\nThe --ref_genome_string argument must be provided. Cannot proceed.\n")
}

### finding reference bed file
ref_bed_path <- list.files(args$ref_dir, pattern = ".*.bed", full.names = TRUE)

# checking only one was found matching pattern search
if ( length(ref_bed_path) != 1 ) { 
    error_message = cat("\nA single reference bed file was not found in the", args$ref_dir, "directory ending in '*bed'. Cannot proceed.\n")
    stop(error_message, call. = FALSE)
}


# getting path to runsheet
runsheet_path <- list.files(args$metadata_dir, pattern = ".*_runsheet.csv", full.names = TRUE)

# checking only one was found matching pattern search
if ( length(runsheet_path) != 1 ) { 
    error_message = cat("\nA single runsheet file was not found in the", args$metadata_dir, "directory ending in '*_runsheet.csv'. Cannot proceed.\n")
    stop(error_message, call. = FALSE)
}


# reading runsheet
runsheet <- read.csv(runsheet_path)

# getting all factors
# adding mock factors so we have more
# runsheet$Factor.Value.Other <- c(rep("Mad", 10), rep("Dog", 6))
# runsheet$Factor.Value.Other2 <- c(rep("LilMac", 6), rep("Fighter", 10))
# runsheet$Factor.Value.Other2 <- c(rep("A", 6), rep("B", 8), rep("C", 2))
factors <- runsheet %>% dplyr::select(starts_with("Factor.Value"))
colnames(factors) = paste("factor", 1:dim(factors)[2], sep = "_")

# checking if there is more than two unique values in a given factor, as methylkit isn't built for that
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
# a lot of this comes from Jonathan Oribello's work in the GeneLab RNAseq workflow
study_df <- factor_df %>% column_to_rownames("sample_id")

# if there are multiple factors, here we are concatenating them to make one combined one
groups <- apply(study_df, 1, paste, collapse = " & ")
group_names <- paste0("(", groups, ")", sep = "") # human readable group names
safe_group_names <- make.names(groups) # group naming compatible with R models
names(safe_group_names) <- group_names

samples_and_combined_factors_df <- data.frame("sample_id" = sample_names, "combined_factor" = safe_group_names)

##### Format contrasts table, defining pairwise comparisons for all groups #####
contrasts <- combn(levels(factor(safe_group_names)), 2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(safe_group_names))), 2)

contrast.names <- c(paste(contrast.names[1,], contrast.names[2,], sep = "v"), 
                    paste(contrast.names[2,], contrast.names[1,], sep = "v"))

contrasts <- cbind(contrasts, contrasts[c(2,1),])
colnames(contrasts) <- contrast.names

if ( args$v ) { cat("\nContrasts table:\n") ; print(contrasts) ; cat("\n") }

# making a single table with info needed for methylkit
sample_meth_info_df <- full_join(samples_and_covs_paths_df, samples_and_combined_factors_df, "sample_id")

# making output directory
dir.create(args$methylkit_output_dir, showWarnings = FALSE)

# reading in transcript features
gene.obj <- readTranscriptFeatures(ref_bed_path, up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, 
                                   unique.prom = TRUE)


####### NOTE TO MIKE DURING DEV #######
# this is about output A vs B AND B vs A, as is currently done with RNAseq
# unlike with deseq2:
    # each contrast needs to be re-run to have the outputs the other way
    # not all the output info can go in one table, so we are doubling all files by reporting both
# so i think ultimately we shouldn't do all contrasts both ways with methylseq data
#######################################

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
            
            curr_treatment_vec <- c(curr_treatment_vec, 0)
            
        } else {
            
            curr_treatment_vec <- c(curr_treatment_vec, 1)
            
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
    
    # getting hyper methylated bases
    # curr_myDiff25p.hyper <- getMethylDiff(curr_myDiff, difference = 25, qvalue = 0.01, type = "hyper")
    
    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_myDiff25p.hyper <- getMethylDiff(curr_myDiff, difference = 1, qvalue = 0.9, type = "hyper")
    ###########################################################
    
    # making table for writing out
    curr_sig_bases_hyper_out_tab <- getData(curr_myDiff25p.hyper) %>% arrange(qvalue)
    
    # getting hypo methylated bases
    # curr_myDiff25p.hypo <- getMethylDiff(curr_myDiff, difference = 25, qvalue = 0.01, type = "hypo")
    
    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_myDiff25p.hypo <- getMethylDiff(curr_myDiff, difference = 1, qvalue = 0.9, type = "hypo")
    ###########################################################
    
    # making table for writing out
    curr_sig_bases_hypo_out_tab <- getData(curr_myDiff25p.hypo) %>% arrange(qvalue)
    
    # getting all differentially methylated bases
    # curr_myDiff25p <- getMethylDiff(curr_myDiff, difference = 25, qvalue = 0.01)
    
    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_myDiff25p <- getMethylDiff(curr_myDiff, difference = 1, qvalue = 0.9)
    ###########################################################
    
    # making table for writing out
    curr_sig_bases_all_out_tab <- getData(curr_myDiff25p) %>% arrange(qvalue)
    
    # writing out tables
    curr_sig_bases_hyper_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-bases.tsv"))
    write.table(curr_sig_bases_hyper_out_tab, curr_sig_bases_hyper_out_tab_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    curr_sig_bases_hypo_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-bases.tsv"))
    write.table(curr_sig_bases_hypo_out_tab, curr_sig_bases_hypo_out_tab_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)

    curr_sig_bases_all_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-bases.tsv"))
    write.table(curr_sig_bases_all_out_tab, curr_sig_bases_all_out_tab_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    ### Tile analysis ###
    # tiling
    curr_tiles_obj <- tileMethylCounts(curr_obj, win.size = 1000, step.size = 1000, cov.bases = 10)
    
    # merging tiled samples
    curr_tiles_meth <- unite(curr_tiles_obj)
    
    # calculating differential methylation on tiles
    curr_tiles_diff <- calculateDiffMeth(curr_tiles_meth, mc.cores = args$mc_cores)
    
    # getting hyper methylated tiles
    # curr_tiles_myDiff25p.hyper <- getMethylDiff(curr_tiles_diff, difference = 25, qvalue = 0.01, type = "hyper")

    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_tiles_myDiff25p.hyper <- getMethylDiff(curr_tiles_diff, difference = 1, qvalue = 0.9, type = "hyper")
    ###########################################################
    
    # making table for writing out
    curr_tiles_sig_hyper_out_tab <- getData(curr_tiles_myDiff25p.hyper) %>% arrange(qvalue)
    
    # getting hypo methylated tiles
    # curr_tiles_myDiff25p.hypo <- getMethylDiff(curr_tiles_diff, difference = 25, qvalue = 0.01, type = "hypo")

    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_tiles_myDiff25p.hypo <- getMethylDiff(curr_tiles_diff, difference = 1, qvalue = 0.9, type = "hypo")
    ###########################################################
    
    # making table for writing out
    curr_tiles_sig_hypo_out_tab <- getData(curr_tiles_myDiff25p.hypo) %>% arrange(qvalue)
    
    # getting all differentially methylated tiles
    # curr_tiles_myDiff25p <- getMethylDiff(curr_tiles_diff, difference = 25, qvalue = 0.01)

    ##### TEST VALUES ONLY - REPLACE WITH ABOVE WHEN DONE #####
    curr_tiles_myDiff25p <- getMethylDiff(curr_tiles_diff, difference = 1, qvalue = 0.9)
    ###########################################################
    
    # making table for writing out
    curr_sig_tiles_all_out_tab <- getData(curr_tiles_myDiff25p) %>% arrange(qvalue)
    
    
    # writing out tables
    curr_tiles_sig_hyper_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypermethylated-tiles.tsv"))
    write.table(curr_tiles_sig_hyper_out_tab, curr_tiles_sig_hyper_out_tab_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)

    curr_tiles_sig_hypo_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-hypomethylated-tiles.tsv"))
    write.table(curr_tiles_sig_hypo_out_tab, curr_tiles_sig_hypo_out_tab_path,
                sep = "\t", quote = FALSE, row.names = FALSE)

    curr_sig_tiles_all_out_tab_path <- file.path(curr_output_dir, paste0(curr_output_prefix, "-sig-diff-methylated-tiles.tsv"))
    write.table(curr_sig_tiles_all_out_tab, curr_sig_tiles_all_out_tab_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    
    
    
    
    
    
    
    
    ### Adding feature information ###
    
    ## adding features to individual-base object
    
    ######## PROBLEM NEXT FEW LINES HERE, MIGHT HAVE TO DO WITH BED FILE NOT BEING SUBSET PROPERLY OR SOMETHING
    curr_diffAnn <- annotateWithGeneParts(as(curr_myDiff25p, "GRanges"), gene.obj)
    
    # making base-level sig table with features 
    curr_sig_all_bases_tab_with_features <- cbind(data.frame(curr_myDiff25p), 
                                                  getAssociationWithTSS(curr_diffAnn), 
                                                  as.data.frame(getMembers(curr_diffAnn))) %>% .[,-c(8)]
    
    write.table(sig_all_bases_tab_with_features, "sig-diff-methylated-bases-with-features.tsv", 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    ## adding features to tiles object
    tiles_diffAnn <- annotateWithGeneParts(as(tiles_myDiff25p, "GRanges"), gene.obj)
    
    # making tiles sig table with features 
    tiles_sig_all_out_tab_with_features <- cbind(data.frame(myDiff25p), 
                                                 getAssociationWithTSS(diffAnn), 
                                                 as.data.frame(getMembers(diffAnn))) %>% .[,-c(8)]
    
    write.table(tiles_sig_all_out_tab_with_features, "sig-diff-methylated-tiles-with-features.tsv", 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    ### Adding functional annotations ###
    # reading in annotation table appropriate for current organism
    # when we have the final location for this information, this will
    # need to be updated to pull a table with the links, rather than
    # the link being hard-coded here in this example
    functional_annots_tab <- 
        read.table("https://figshare.com/ndownloader/files/35939642", sep = "\t", 
                   quote = "", header = TRUE)
    
    # reading in gene to transcript mapping file
    gene_transcript_map <- 
        read.table("subset-test-results/Mus_musculus.GRCm39.107-gene-to-transcript-map.tsv", sep = "\t", 
                   col.names = c("gene_ID", "feature.name"))
    
    
    ## for individual-base output
    # for each transcript ID in the sig_all_bases_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table
    sig_all_bases_tab_with_features_and_gene_IDs <- 
        left_join(sig_all_bases_tab_with_features, gene_transcript_map)
    
    # now adding full annotations
    sig_all_bases_tab_with_features_and_annots <- 
        left_join(sig_all_bases_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
    
    # and writing out
    write.table(sig_all_bases_tab_with_features_and_annots, 
                "sig-diff-methylated-bases-with-features-and-annots.tsv", 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    # making and writing out a table of base-level percent methylated
    perc.meth <- percMethylation(meth, rowids = TRUE)
    write.table(perc.meth, "base-level-percent-methylated.tsv", sep = "\t", 
                quote = FALSE, row.names = TRUE, col.names=NA)
    
    
    ## for tiles output
    # for each transcript ID in the tiles_sig_all_out_tab_with_features table, getting 
    # its corresponding gene ID and adding that to the table
    sig_all_tiles_tab_with_features_and_gene_IDs <- 
        left_join(tiles_sig_all_out_tab_with_features, gene_transcript_map)
    
    # now adding full annotations
    sig_all_tiles_tab_with_features_and_annots <- 
        left_join(sig_all_tiles_tab_with_features_and_gene_IDs, 
                  functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
    
    # and writing out
    write.table(sig_all_tiles_tab_with_features_and_annots, 
                "sig-diff-methylated-tiles-with-features-and-annots.tsv", 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    # making and writing out a table of tile-level percent methylated
    tiles_perc.meth <- percMethylation(tiles_meth, rowids = TRUE)
    write.table(tiles_perc.meth, "tile-level-percent-methylated.tsv", sep = "\t", 
                quote = FALSE, row.names = TRUE, col.names=NA)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}


sample_names <- c("M21", "M22", "M23")
bismark_cov_paths

mod_paths_order <- bismark_cov_paths[c(3,1,2)]



# # 
# # ## beginning methylkit
# # # reading into memory
# # myobj <- methRead(location = file.list,
# #                   sample.id = sample.list,
# #                   assembly = "Mmus_GRCm39",
# #                   pipeline = "bismarkCoverage",
# #                   header = FALSE,
# #                   treatment = c(1,1,1,0,0,0),
# #                   mincov = 10)
# 
# 


### helper functions ###
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

