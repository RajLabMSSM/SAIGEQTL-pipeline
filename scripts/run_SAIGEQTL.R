options(echo=TRUE)
## RUN SAIGEQTL

# Beomjin Jang
# 2026
# wraps around the SAIGEQTL binary
# for massive parallelisation

## Arguments:
# all the paths to the input file lists needed
# path to phenotype metadata
# chunk number (n) and chunk number (i of n)

library(parallel)
library(tidyverse)
library(dplyr)
library(optparse)

# Define and parse command-line options
option_list <- list(
  make_option(c('--chrom'), help = "which chromosome to run on"),
  make_option(c('--threads'), help = "number of threads for parallelization", default = 4L),
  make_option(c('--pheno_meta'), help = 'phenotype metadata file'),
  make_option(c('--chunk_total', '-n'), help = 'total number of chunks', default = 1L),
  make_option(c('--chunk', '-i'), help = 'current chunk number', default = 1L),
  make_option(c('--prefix'), help = 'stem of out file path'),
  make_option(c('--geno_file'), help = "path to PLINK file prefix" ),
  make_option(c('--pheno_file'), help = "path to phenotype file"),
  make_option(c('--cov_cols'), default = '', help = "comma separated list of covariates"),
  make_option(c('--sample_cov_cols'), default = '', help = "comma separated list of sample-level covariates"),
  make_option(c('--sample_cols'), help = "column name for sample IDs in phenotype file"),
  make_option('--offset_col', type='character', default=''),
  make_option(c('--SAIGEQTL'), help = "full path to SAIGEQTL executable"),
  make_option('--trait_type', default = 'count'),
  make_option('--qcov_cols', default = ''),
  make_option('--tol', type = 'double', default = 0.02),
  make_option('--cis_window', type = 'double', default = 1e6),
  make_option('--minMAF', type = 'double', default = .01),
  make_option('--minMAC', type = 'double', default = 20),
  make_option('--SPAcutoff', type = 'double', default = 2),
  make_option('--markers_per_chunk', type = 'integer', default = 10000L)
)

# option.parser <- OptionParser(option_list=option_list)
# opt <- parse_args(option.parser)
opt <- parse_args(OptionParser(option_list = option_list))
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

num_cores <- opt$threads
message("Using ", num_cores, " cores for parallel processing")

chrom <- opt$chrom
meta_file <- absPath(opt$pheno_meta)
geno_file <- sub('[.]bed$', '', absPath(paste0(opt$geno_file, '.bed')))
pheno_file <- absPath(opt$pheno_file)
i_chunk <- as.numeric(opt$chunk)
n_chunk <- as.numeric(opt$chunk_total)
prefix <- absPath(opt$prefix)
saigeqtl_bin <- opt$SAIGEQTL

# Capture covariate lists from arguments
cov_cols <- opt$cov_cols
sample_cov_cols <- opt$sample_cov_cols
sample_cols <- opt$sample_cols
cis_window <- opt$cis_window

# --- Directory Setup ---
# if (!dir.exists(prefix)) {
#   dir.create(prefix, recursive = TRUE)
# }
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

# --- Function Definitions ---

# Splits the list of features into a given number of chunks (n)
split_chunks <- function(meta, i, n){
  if (n < 1 || i < 1 || i > n || n > nrow(meta)) stop('Invalid chunk count')
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
  
  chromosome <- meta_j$chr
  
  regionFile <- paste0(prefix_region, '/cis_region_', feature_j, '.txt')
  write_tsv(meta_j[,c('chr','start','end')], regionFile, col_names = FALSE)
  
  prefix_step_01_feature_j <- paste0(prefix_step_01, chromosome,'_',feature_j)
  prefix_step_02_feature_j <- paste0(prefix_step_02, '/SAIGEQTL_', chromosome, '_', feature_j, '.txt')
  # Skip this gene if its Step 2 result already exists.
  if (file.exists(prefix_step_02_feature_j)) {
    message("Skipping feature: ", feature_j, " (result exists)")
    return(TRUE)
  }
  
  # prefix_step_03_feature_j <- paste0(prefix_step_03, '/SAIGEQTL_',chromosome,'_', feature_j,'txt')
  
  # Define paths to Step 1 output files
  rda_file <- paste0(prefix_step_01_feature_j, ".rda")
  var_ratio_file <- paste0(prefix_step_01_feature_j, ".varianceRatio.txt")
  
  geno_file_step1 <- paste0(geno_file, '_step1')
  
  # Command for Step 1: Fit Null GLMM
  cmd1 <- paste0( "unset R_LIBS_USER; unset R_LIBS_SITE; unset XDG_RUNTIME_DIR; ml singularity; singularity exec ", shQuote(saigeqtl_bin), " step1_fitNULLGLMM_qtl.R",
                  " --useSparseGRMtoFitNULL=FALSE",
                  " --useGRMtoFitNULL=FALSE",
                  " --phenoFile=", shQuote(pheno_file),
                  " --phenoCol=", shQuote(feature_j),
                  " --covarColList=", shQuote(cov_cols),
                  " --sampleCovarColList=", shQuote(sample_cov_cols),
                  " --sampleIDColinphenoFile=", shQuote(sample_cols),
                  " --traitType=", opt$trait_type,
                  " --outputPrefix=", shQuote(prefix_step_01_feature_j),
                  " --skipVarianceRatioEstimation=FALSE",
                  " --isRemoveZerosinPheno=FALSE",
                  " --isCovariateOffset=FALSE",
                  " --isCovariateTransform=TRUE",
                  " --skipModelFitting=FALSE",
                  " --tol=", opt$tol,
                  " --plinkFile=", shQuote(geno_file_step1),
                  " --IsOverwriteVarianceRatioFile=TRUE"
  )
  if (nzchar(opt$offset_col)) {
    cmd1 <- paste0(cmd1, " --offsetCol=", shQuote(opt$offset_col))
  }
  
  if (nzchar(opt$qcov_cols)) {
    cmd1 <- paste0(cmd1, " --qCovarColList=", shQuote(opt$qcov_cols))
  }
  
  # Execute Step 1 and capture its exit status
  message("Running Step 1 for feature: ", feature_j)
  print(cmd1)
  cmd1_status <- system(cmd1)
  
  # Check if Step 1 was successful AND created the required files
  if (cmd1_status == 0 && file.exists(rda_file) && file.exists(var_ratio_file)) {
    message("Step 1 successful for ", feature_j, ". Proceeding to Step 2.")
    
    # Command for Step 2: Main association testing
    cmd2 <- paste0( "unset R_LIBS_USER; unset R_LIBS_SITE; unset XDG_RUNTIME_DIR; ml singularity; singularity exec ", shQuote(saigeqtl_bin), " step2_tests_qtl.R",
                    " --bedFile=", shQuote(paste0(geno_file,".bed")),
                    " --bimFile=", shQuote(paste0(geno_file,".bim")),
                    " --famFile=", shQuote(paste0(geno_file,".fam")),
                    " --SAIGEOutputFile=", shQuote(prefix_step_02_feature_j),
                    " --chrom=", chromosome,
                    " --minMAF=", opt$minMAF,
                    " --minMAC=", opt$minMAC,
                    " --LOCO=FALSE",
                    " --GMMATmodelFile=", shQuote(rda_file),
                    " --SPAcutoff=", opt$SPAcutoff,
                    " --varianceRatioFile=", shQuote(var_ratio_file),
                    " --rangestoIncludeFile=", shQuote(regionFile),
                    " --markers_per_chunk=", opt$markers_per_chunk
    )
    
    # Execute Step 2
    print(cmd2)
    system(cmd2)
    
  } else {
    # If Step 1 failed or didn't produce output, print a warning and stop.
    warning("Step 1 failed for feature: ", feature_j, ". Exit status: ", cmd1_status, ". Skipping Step 2.")
  }

}

# --- Main Execution ---

# Read and prepare metadata
meta <- readr::read_tsv(meta_file, col_names = c("chr", "start", "end", "feature"),
                        col_types = "cddc")

# Define the final output file expected by Snakemake
final_out_file <- paste0(prefix, "/SAIGEQTL_", chrom, "_chunk_", i_chunk)
if (file.exists(final_out_file)) unlink(final_out_file)
if (!(chrom %in% meta$chr)) {
  stop("Chromosome ", chrom, " not in meta file")
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
if (!file.create(final_out_file)) stop('Could not create completion file')
