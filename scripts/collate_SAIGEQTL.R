## Beomjin Jang, Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)

option_list <- list(
   make_option(c('--output_file', '-o'), help = 'name of out file', default = "results/example/example"),
   make_option(c('--chrom'), help = 'the chromosome', default = ""),
   make_option(c('--prefix', '-p'), help = 'prefix of chr specific top assoc files', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix
chrom <- opt$options$chrom

# output_file = "/sc/arion/projects/bigbrain/SingleBrain/SAIGEQTL-pipeline/results/MG/chr22_full_assoc.tsv.gz"
# prefix = '/sc/arion/projects/bigbrain/SingleBrain/SAIGEQTL-pipeline/results/MG/SAIGEQTL_tmp/'
# chrom='chr22'
prefix_step_02 <- paste0(prefix, '/step_02/')

# top_pattern <- paste0("SAIGEQTL_^chr[0-9]+_ENSG_")
top_pattern <- paste0("SAIGEQTL_",chrom,'_')

print(prefix_step_02)
inputs <- list.files(prefix_step_02, pattern = top_pattern, recursive = FALSE, full.names = TRUE)
inputs <- gsub('.index','',inputs) %>% unique()

message(" * Processing full associations for QTL peak ")

# Check for input files
if (length(inputs) == 0) {
  message(" * WARNING: No top files found for peak Writing an empty file.")
  write_tsv(tibble(), output_file)  # Write an empty output file and exit
  quit(save = "no")
}

# Read and combine input files
process_sumstats <- function(x){
   feature_j <-str_extract(x, "ENSG\\d+\\.\\d+")

   sumstats_feature_j <- read_tsv(x)
   if(nrow(sumstats_feature_j)==0) { return(NULL) } else{
      message(feature_j)
      sumstats_feature_j$feature <- feature_j
      
      # Columns to move to the front
      cols_to_front <- c("feature")
      
      # Get the names of the remaining columns
      remaining_cols <- setdiff(names(sumstats_feature_j), cols_to_front)
      
      # Combine them in the desired order
      new_order <- c(cols_to_front, remaining_cols)
      
      # Reorder the data frame
      sumstats_feature_j <- sumstats_feature_j[, new_order]
      
      sumstats_feature_j <- sumstats_feature_j[sumstats_feature_j$p.value != 0, ]
      sumstats_feature_j$p_bonf <- p.adjust(sumstats_feature_j$p.value, method = "bonferroni")
      sumstats_feature_j$p_FDR <- p.adjust(sumstats_feature_j$p.value, method = "BH" )
      
      return(sumstats_feature_j)
   }

}

res <- purrr::map_df(inputs, ~{
   process_sumstats(.x)
}, ) 

output_file_unzip <- gsub('_full_assoc.tsv.gz' ,'_full_assoc.tsv', output_file)
header_file <- paste0(output_file_unzip, "_header")

if( chrom == "chr1" ){
   # write column names as separate file
   write_tsv(res[0,], path = header_file, col_names = TRUE)
   write_tsv(res, path = output_file_unzip, col_names = FALSE) 
}else{
   write_tsv(res, path = output_file_unzip, col_names = FALSE)
}

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
