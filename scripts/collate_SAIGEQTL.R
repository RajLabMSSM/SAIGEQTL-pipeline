## Beomjin Jang, Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)
library(data.table)

option_list <- list(
  make_option(c('--output_file', '-o'), help = 'name of out file', default = "results/example/example"),
  make_option(c('--chrom'), help = 'the chromosome', default = ""),
  make_option(c('--prefix', '-p'), help = 'prefix of chr specific top assoc files', default = "results/example/example"),
  make_option(c('--pheno_meta'), help = 'four-column phenotype metadata without header')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix
chrom <- opt$options$chrom
pheno_meta <- opt$options$pheno_meta

prefix_step_02 <- paste0(prefix, '/step_02/')

meta <- fread(pheno_meta, header=FALSE,
              col.names=c('chr','start','end','feature'), colClasses=c(V1='character'))
features <- meta[chr == chrom, unique(feature)]
inputs <- paste0(prefix_step_02, 'SAIGEQTL_', chrom, '_', features, '.txt')
if (!length(inputs)) stop('No features for chromosome ', chrom)
if (any(!file.exists(inputs))) stop('Missing Step 2 files: ', paste(inputs[!file.exists(inputs)], collapse=', '))

message(" * Processing full associations for QTL peak ")
# Read and combine input files
process_sumstats <- function(x, feature_j){
  sumstats_feature_j <- readr::read_tsv(
    x,
    col_types = readr::cols(
      CHR           = readr::col_character(),
      POS           = readr::col_double(),
      MarkerID      = readr::col_character(),
      Allele1       = readr::col_character(),
      Allele2       = readr::col_character(),
      AC_Allele2    = readr::col_double(),
      AF_Allele2    = readr::col_double(),
      MissingRate   = readr::col_double(),
      BETA          = readr::col_double(),
      SE            = readr::col_double(),
      Tstat         = readr::col_double(),
      var           = readr::col_double(),
      p.value       = readr::col_double(),
      p.value.NA    = readr::col_double(),
      Is.SPA        = readr::col_logical(),
      N             = readr::col_double()
    )
  )
  if (!all(c('CHR','POS','MarkerID','p.value') %in% names(sumstats_feature_j))) stop('Unexpected schema: ', x)
  sumstats_feature_j <- sumstats_feature_j %>%
    dplyr::mutate(
      feature = feature_j,
      p_bonf = p.adjust(p.value, method = "bonferroni"),
      p_FDR = p.adjust(p.value, method = "BH")
    ) %>%
    dplyr::relocate(feature)
  
  return(sumstats_feature_j)
}

res <- purrr::map2_df(inputs, features, function(x, feature_j) {
  process_sumstats(x, feature_j)
})

output_file_unzip <- gsub('_full_assoc.tsv.gz' , '_full_assoc.tsv', output_file)
header_file <- paste0(output_file_unzip, "_header")

# write column names as separate file
write_tsv(res[0, ], path = header_file, col_names = TRUE)
write_tsv(res, path = output_file_unzip, col_names = FALSE)

# bgzip
message( "* bgzipping " )
bgzip_cmd <- paste0(" ml tabix; bgzip -f -c ", output_file_unzip, " > ", output_file)
message( " * ", bgzip_cmd)
system(bgzip_cmd)

# Removing
message( "* Removing " )
remove_cmd <- paste0(" rm ", output_file_unzip)
message( " * ", remove_cmd)
system(remove_cmd)

message(" * Full associations for peak written to: ", output_file)
message(" * Collation complete.")
