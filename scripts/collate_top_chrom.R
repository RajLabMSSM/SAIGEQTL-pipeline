## Beomjin Jang, Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)
library(data.table)
library(readr)

option_list <- list(
  make_option(c('--output_file', '-o'), help = 'name of out file'),
  make_option(c('--prefix', '-p'), help = 'CASTIE temporary folder'),
  make_option(c('--step3'), help = 'combined Step 3 long-format file'),
  make_option(c('--step4_prefix'), help = 'Step 4 path and filename prefix, e.g. results/Mic/step4_Mic_')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix
step3 <- opt$options$step3
step4_prefix <- opt$options$step4_prefix

# output_file = "/sc/arion/projects/bigbrain/SingleBrain/SAIGEQTL-pipeline/results/MG/MG_top_assoc.tsv.gz"
# prefix = '/sc/arion/projects/bigbrain/SingleBrain/CASTIE-pipeline/results/Mic_test/CASTIE_tmp/'
# step3 = '/sc/arion/projects/bigbrain/SingleBrain/CASTIE-pipeline/results/Mic_test/step3_longformat.txt'
prefix_step_02 <- paste0(prefix, '/step_02/')


# Read all tests/contexts, then calculate SNP corrections once per gene.
res_top <- as.data.table(readr::read_tsv(step3, show_col_types=FALSE))
if (!nrow(res_top)) stop('No top associations in Step 3')
features <- unique(res_top$Gene)
inputs <- list.files(prefix_step_02, pattern='^SAIGEQTL_.*[.]txt$', full.names=TRUE)

extract_sumstats <- function(j){
  feature_j <- features[j]
  input_j <- inputs[sub('^SAIGEQTL_[^_]+_', '', sub('[.]txt$', '', basename(inputs))) == feature_j]
  if (length(input_j) != 1) stop('Expected one Step 2 file for ', feature_j)
  sumstats_feature_j <- as.data.table(readr::read_tsv(
    input_j,
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
    ),
    show_col_types = FALSE
  ))
  sumstats_feature_j[, p_bonf := p.adjust(p.value, method='bonferroni')]
  sumstats_feature_j[, p_FDR := p.adjust(p.value, method='BH')]
  
  # A SNP may be the lead for several contexts; keep one row for each test.
  top_j <- res_top[Gene == feature_j]
  if (anyDuplicated(sumstats_feature_j$MarkerID)) stop('Duplicate SNP IDs for ', feature_j)
  index <- match(top_j$top_MarkerID, sumstats_feature_j$MarkerID)
  if (anyNA(index)) stop('Missing lead SNP for ', feature_j)
  sumstats_feature_j <- sumstats_feature_j[index]
  sumstats_feature_j[, feature := feature_j]
  sumstats_feature_j[, ACAT_p := top_j$ACAT_p]
  sumstats_feature_j[, pval_column := top_j$pval_column]
  setcolorder(sumstats_feature_j, c('feature', setdiff(names(sumstats_feature_j), 'feature')))
  return(sumstats_feature_j)
}
res <- purrr::map_df(seq_along(features), extract_sumstats)
res <- dplyr::arrange(res, ACAT_p)

# Attach Step 4 q-values
egene_file <- paste0(step4_prefix, 'pval_main_egene.tsv')

egene_qv <- readr::read_tsv(
  egene_file,
  show_col_types = FALSE,
  col_types = cols(
    .default = col_skip(),
    Gene = col_character(),
    pval_column = col_character(),
    top_MarkerID = col_character(),
    qv = col_double()
  )
)

res <- dplyr::left_join(
  res,
  egene_qv,
  by = c(
    'feature' = 'Gene',
    'pval_column' = 'pval_column',
    'MarkerID' = 'top_MarkerID'
  )
)
res <- dplyr::relocate(res, pval_column, .after = dplyr::last_col())

# Save main associations
output_file <- paste0(output_file, '_pval_main.tsv.gz')
output_file_unzip <- gsub('[.]tsv[.]gz$', '.tsv', output_file)

write_tsv(res, output_file_unzip)

# bgzip
message("* bgzipping ")
bgzip_cmd <- paste0("bgzip -f -c ",
                    shQuote(output_file_unzip),
                    " > ",
                    shQuote(output_file))
message(" * ", bgzip_cmd)
if (system(bgzip_cmd) != 0)
  stop('bgzip failed')

# Removing
message("* Removing ")
remove_cmd <- paste0("rm ", shQuote(output_file_unzip))
message(" * ", remove_cmd)
system(remove_cmd)

message(" * Top associations for peak written to: ", output_file)

