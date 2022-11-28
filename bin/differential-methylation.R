library(optparse)
library(tidyverse)

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


# directory of metadata
parser <- add_option(parser, c("--limit_samples_to"),
                     default = "all",
                     help = "Limits the number of samples being processed (won't do real factor comparisons)")


args <- parse_args(parser)


# getting path to runsheet
runsheet_filename <- list.files(args$metadata_dir, pattern = ".*_runsheet.csv")

# checking only one was found matching pattern search
if ( length(runsheet_filename) != 1 ) { 
    error_message = cat("More than one file in ", args$metadata_dir, " found ending in '*_runsheet.csv'. Cannot proceed.\n")
    stop(error_message, call. = FALSE)
}

# reading runsheet
runsheet_path <- file.path(args$metadata, runsheet_filename)

runsheet <- read.csv(runsheet_path)

# getting all factors
# adding mock factor so we have 2
# runsheet$Factor.Value.Other <- c(rep("Mike", 8), rep("Lee", 8))
factors <- runsheet %>% select(starts_with("Factor.Value"))
colnames(factors) = paste("factor", 1:dim(factors)[2], sep = "_")

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

if ( args$v ) { cat("\nFactor df:\n") ; factor_df }

# getting list of *.bismark.cov.gz files
bismark_cov_filenames <- list.files(args$bismark_methylation_calls_dir, pattern = ".*.bismark.cov.gz")

bismark_cov_paths <- file.path(args$bismark_methylation_calls_dir, bismark_cov_filenames)

# making sure files list matches length of runsheet
if ( dim(runsheet)[1] != length(bismark_cov_paths) ) {

    error_message = cat("The number of '*.bismark.cov.gz' files found at ", args$bismark_methylation_calls_dir, " does not match the number of samples specified in the runsheet. Cannot proceed.\n")
    stop(error_message, call. = FALSE)
    
}

if ( args$v ) { cat("\nBismark coverage files:\n") ; bismark_cov_paths }


## reading ONE methylseq object into memory looks like this
myobj <- methRead(location = file.list,
                  sample.id = sample.list,
                  assembly = "Mmus_GRCm39",
                  pipeline = "bismarkCoverage",
                  header = FALSE,
                  treatment = c(1,1,1,0,0,0),
                  mincov = 10)

# we want to loop for the number of things we have


factor_df
study_df <- factor_df %>% column_to_rownames("sample_id")

##### Format groups and indicate the group that each sample belongs to #####
# concatenate multiple factors into one condition per sample
groups <- apply(study_df, 1, paste, collapse = " & ")

group_names <- paste0("(", groups, ")", sep = "") # human readable group names
safe_group_names <- make.names(groups) # group naming compatible with R models
names(safe_group_names) <- group_names

##### Format contrasts table, defining pairwise comparisons for all groups #####
contrasts <- combn(levels(factor(safe_group_names)), 2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(safe_group_names))), 2)

contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names)


## Create data frame defining which group each sample belongs to
sampleTable <- data.frame(condition = factor(safe_group_names))
row.names(sampleTable) <- factor_df$sample_id
factor(sampleTable$condition)



## beginning methylkit
# reading into memory
myobj <- methRead(location = file.list,
                  sample.id = sample.list,
                  assembly = "Mmus_GRCm39",
                  pipeline = "bismarkCoverage",
                  header = FALSE,
                  treatment = c(1,1,1,0,0,0),
                  mincov = 10)


