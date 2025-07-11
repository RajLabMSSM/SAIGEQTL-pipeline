## Beomjin Jang, Jack Humphrey, Kailash BP 2025
## Collate top associations per peak

library(optparse)
library(tidyverse)

option_list <- list(
   make_option(c('--output_file', '-o'), help = 'name of out file', default = "results/example/example"),
   make_option(c('--prefix', '-p'), help = 'prefix of chr specific top assoc files', default = "results/example/example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
output_file <- opt$options$output_file
prefix <- opt$options$prefix

# output_file = "/sc/arion/projects/bigbrain/SingleBrain/SAIGEQTL-pipeline/results/MG/MG_top_assoc.tsv.gz"
# prefix = '/sc/arion/projects/bigbrain/SingleBrain/SAIGEQTL-pipeline/results/MG/SAIGEQTL_tmp/'
prefix_step_02 <- paste0(prefix, '/step_02/')
prefix_step_03 <- paste0(prefix, '/step_03/')


# top_pattern <- paste0("SAIGEQTL_^chr[0-9]+_ENSG_")
top_pattern <- paste0("SAIGEQTL_")

print(prefix_step_03)
inputs <- list.files(prefix_step_03, pattern = top_pattern, recursive = FALSE, full.names = TRUE)

message(" * Processing top associations for QTL peak ")

# Check for input files
if (length(inputs) == 0) {
  message(" * WARNING: No top files found for peak Writing an empty file.")
  write_tsv(tibble(), output_file)  # Write an empty output file and exit
  quit(save = "no")
}

# Read and combine input files
res_top <- map_df(inputs, read_tsv)
colnames(res_top) <- c('feature','ACAT_p','MarkerID','p.value')

extract_sumstats <- function(x, df_top){
   feature_j <-str_extract(x, "ENSG\\d+\\.\\d+")
   message(feature_j)
   df_top_feautre_j <-  df_top[df_top$feature == feature_j,]
   
   path_sumstats <- sumstats_feature_j <- gsub(prefix_step_03, prefix_step_02, x )
   sumstats_feature_j <- read_tsv(path_sumstats)
   
   sumstats_feature_j$p_bonf <- p.adjust(sumstats_feature_j$p.value, method = "bonferroni")
   sumstats_feature_j$p_FDR <- p.adjust(sumstats_feature_j$p.value, method = "BH" )
   
   sumstats_feature_j <- sumstats_feature_j[sumstats_feature_j$MarkerID %in% df_top_feautre_j$MarkerID, ]
   sumstats_feature_j$feature <- feature_j
   
   return(sumstats_feature_j)
}
res <- purrr::map_df(inputs, ~{
   extract_sumstats(.x, res_top)
}, ) 

res <- right_join(res, res_top, by=c('feature','MarkerID','p.value'))
res <- res %>% dplyr::filter(ACAT_p!=0) %>% 
   arrange(ACAT_p) 


# Columns to move to the front
cols_to_front <- c("feature")

# Get the names of the remaining columns
remaining_cols <- setdiff(names(res), cols_to_front)

# Combine them in the desired order
new_order <- c(cols_to_front, remaining_cols)

# Reorder the data frame
res <- res[, new_order]
print(res)

output_file_unzip <- gsub('_top_assoc.tsv.gz' ,'_top_assoc.tsv', output_file)

write_tsv(res, output_file_unzip)

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

message(" * Top associations for peak written to: ", output_file)

