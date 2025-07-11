options(echo=TRUE)
## RUN SAIGEQTL

# Beomjin Jang & Kailash BP & Jack Humphrey
# 2025
# wraps around the SAIGEQTL binary
# for massive parallelisation

## Arguments:
# all the paths to the input file lists needed
# path to phenotype metadata
# chunk number (n) and chunk number (i of n)

library(parallel)
library(tidyverse)
library(optparse)

# Define and parse command-line options
option_list <- list(
  make_option(c('--chrom'), help = "which chromosome to run on"),
  make_option(c('--threads'), help = "number of threads for parallelization", default = 4),
  make_option(c('--pheno_meta'), help = 'phenotype metadata file'),
  make_option(c('--chunk_total', '-n'), help = 'total number of chunks', default = 10),
  make_option(c('--chunk', '-i'), help = 'current chunk number', default = 1),
  make_option(c('--prefix'), help = 'stem of out file path'),
  make_option(c('--geno_file'), help = "path to PLINK file prefix" ),
  make_option(c('--pheno_file'), help = "path to phenotype file"),
  make_option(c('--cov_cols'), help = "comma separated list of covariates"),
  make_option(c('--sample_cov_cols'), help = "comma separated list of sample-level covariates"),
  make_option(c('--sample_cols'), help = "column name for sample IDs in phenotype file"),
  make_option(c('--SAIGEQTL'), help = "full path to SAIGEQTL executable")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

# --- Argument Processing ---

# Function to get the absolute path of a file
absPath <- function(path){
  if(is.null(path)){ return(NULL) }
  # Check if the path is a file or a directory that exists
  if(file.exists(path) || dir.exists(path)){
    return(normalizePath(path))
  }
  return(path)
}

num_cores <- as.numeric(opt$threads)
message("Using ", num_cores, " cores for parallel processing")

chrom <- opt$chrom
meta_file <- absPath(opt$pheno_meta)
geno_file <- absPath(opt$geno_file)
pheno_file <- absPath(opt$pheno_file)
i_chunk <- as.numeric(opt$chunk)
n_chunk <- as.numeric(opt$chunk_total)
prefix <- absPath(opt$prefix)
saigeqtl_bin <- opt$SAIGEQTL

# Capture covariate lists from arguments
cov_cols <- opt$cov_cols
sample_cov_cols <- opt$sample_cov_cols
sample_cols <- opt$sample_cols

cis_window <- 1e6

# --- Directory Setup ---
if (!dir.exists(prefix)) {
  dir.create(prefix, recursive = TRUE)
}
prefix_step_01 <- paste0(prefix, '/step_01/')
if (!dir.exists(prefix_step_01)) {
  dir.create(prefix_step_01, recursive = TRUE)
}
prefix_region <- paste0(prefix, '/cis_region/')
if (!dir.exists(prefix_region)) {
  dir.create(prefix_region, recursive = TRUE)
}
prefix_step_02 <- paste0(prefix, '/step_02/')
if (!dir.exists(prefix_step_02)) {
  dir.create(prefix_step_02, recursive = TRUE)
}
prefix_step_03 <- paste0(prefix, '/step_03/')
if (!dir.exists(prefix_step_03)) {
  dir.create(prefix_step_03, recursive = TRUE)
}

# --- Function Definitions ---

# Splits the list of features into a given number of chunks (n)
split_chunks <- function(meta, i, n){
  bins <- dplyr::ntile(1:nrow(meta), n)
  chunk <- meta[bins == i, ]
  return(chunk)
}

# Main function to construct and run SAIGEQTL commands
run_SAIGEQTL <- function(meta_loc, j){
  meta_j <- meta_loc[j, ]
  feature_j <- meta_j$feature
  
  meta_j$start <- as.numeric(meta_j$start) - cis_window
  if(is.na(meta_j$start) || meta_j$start < 0){ meta_j$start <- 0 }
  meta_j$end <- as.numeric(meta_j$end) + cis_window
  
  chromosome <- as.numeric(gsub('chr','',meta_j$chr))
  chromosome <- meta_j$chr
  
  regionFile <- paste0(prefix_region, '/cis_region_', feature_j, '.txt')
  write_tsv(meta_j[,c('chr','start','end')], regionFile, col_names = FALSE)
  
  prefix_step_01_feature_j <- paste0(prefix_step_01, chromosome,'_',feature_j)
  prefix_step_02_feature_j <- paste0(prefix_step_02, '/SAIGEQTL_',chromosome,'_', feature_j)
  prefix_step_03_feature_j <- paste0(prefix_step_03, '/SAIGEQTL_',chromosome,'_', feature_j)
  
  # Define paths to Step 1 output files
  rda_file <- paste0(prefix_step_01_feature_j, ".rda")
  var_ratio_file <- paste0(prefix_step_01_feature_j, ".varianceRatio.txt")

  geno_file_step1 <- paste0(geno_file, '_step1')
    
  # Command for Step 1: Fit Null GLMM
  cmd1 <- paste0( "unset R_LIBS_USER; unset R_LIBS_SITE; unset XDG_RUNTIME_DIR; ml singularity; singularity exec ", saigeqtl_bin, " step1_fitNULLGLMM_qtl.R",
                  " --useSparseGRMtoFitNULL=FALSE",
                  " --useGRMtoFitNULL=FALSE",
                  " --phenoFile=", pheno_file,
                  " --phenoCol=", feature_j,
                  " --covarColList=", cov_cols,
                  " --sampleCovarColList=", sample_cov_cols,
                  " --sampleIDColinphenoFile=", sample_cols,
                  " --traitType=count",
                  " --outputPrefix=", prefix_step_01_feature_j,
                  " --plinkFile=", geno_file_step1,
                  " --IsOverwriteVarianceRatioFile=TRUE"
  )
  
  # Execute Step 1 and capture its exit status
  message("Running Step 1 for feature: ", feature_j)
  print(cmd1)
  cmd1_status <- system(cmd1)
  
  # Check if Step 1 was successful AND created the required files
  if (cmd1_status == 0 && file.exists(rda_file) && file.exists(var_ratio_file)) {
    message("Step 1 successful for ", feature_j, ". Proceeding to Step 2.")

    # Command for Step 2: Main association testing
    cmd2 <- paste0( "unset R_LIBS_USER; unset R_LIBS_SITE; unset XDG_RUNTIME_DIR; ml singularity; singularity exec ", saigeqtl_bin, " step2_tests_qtl.R",
                    " --bedFile=", geno_file,".bed",
                    " --bimFile=", geno_file,".bim",
                    " --famFile=", geno_file,".fam",
                    " --SAIGEOutputFile=", prefix_step_02_feature_j,
                    " --chrom=", chromosome,
                    " --minMAF=0.01",
                    " --minMAC=20",
                    " --LOCO=FALSE",
                    " --GMMATmodelFile=", rda_file,
                    " --SPAcutoff=2",
                    " --varianceRatioFile=", var_ratio_file,
                    " --rangestoIncludeFile=", regionFile,
                    " --markers_per_chunk=10000"
    )

    # Execute Step 2
    print(cmd2)
    system(cmd2)

  } else {
    # If Step 1 failed or didn't produce output, print a warning and stop.
    warning("Step 2 failed for feature: ", feature_j, ". Exit status: ", cmd1_status, ". Skipping Step 2.")
  }
  
  # Check if Step 2 was successful AND created the required files
  if (file.exists(prefix_step_02_feature_j)) {
    message("Step 2 successful for ", feature_j, ". Proceeding to Step 3.")
    
    # Command for Step 2: Main association testing
    cmd3 <- paste0( "unset R_LIBS_USER; unset R_LIBS_SITE; unset XDG_RUNTIME_DIR; ml singularity; singularity exec ", saigeqtl_bin, " step3_gene_pvalue_qtl.R",
                    " --assocFile=", prefix_step_02_feature_j,
                    " --geneName=",feature_j,
                    " --genePval_outputFile=", prefix_step_03_feature_j
    )
    
    # Execute Step 2
    print(cmd3)
    system(cmd3)
    
  } else {
    # If Step 1 failed or didn't produce output, print a warning and stop.
    warning("Step 3 failed for feature: ", feature_j, ". Exit status: Skipping Step 3.")
  }
}

# --- Main Execution ---

# Read and prepare metadata
meta <- readr::read_tsv(meta_file, col_names = c("chr", "start", "end", "feature"))

# Define the final output file expected by Snakemake
final_out_file <- paste0(prefix, "/SAIGEQTL_", chrom, "_chunk_", i_chunk)

if (!(chrom %in% meta$chr)) {
  warning("WARNING: Chromosome ", chrom, " not in meta file.")
  file.create(final_out_file) # Create empty sentinel file for Snakemake
  quit(save = "no")
}

# Filter metadata for the specified chromosome and chunk
meta_loc <- meta[meta$chr == chrom,]
meta_loc <- split_chunks(meta_loc, i_chunk, n_chunk)

message(" * ", nrow(meta_loc), " features in chunk ", i_chunk, " of ", n_chunk)

# Run SAIGEQTL in parallel if features exist in the chunk
if (nrow(meta_loc) > 0) {
  mclapply(1:nrow(meta_loc), function(j) {
    run_SAIGEQTL(meta_loc, j)
  }, mc.cores = num_cores)
}

# Create the final output file to signal completion to Snakemake
file.create(final_out_file)