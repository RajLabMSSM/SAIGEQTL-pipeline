#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qvalue)
  library(optparse)
})

option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to Step 3 long-format input file [required]"),
  make_option("--outdir", type = "character", default = ".",
              help = "Output directory [default: current directory]"),
  make_option("--prefix", type = "character", default = "",
              help = "Optional output file prefix [default: none]"),
  make_option("--fdr", type = "numeric", default = 0.05,
              help = "FDR threshold for eGene calling [default: 0.05]"),
  make_option("--main_col", type = "character", default = "pval_main",
              help = "Column name for main eQTL p-values [default: pval_main]")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input)) stop("--input is required")
if (!file.exists(opt$input)) stop("Input file does not exist: ", opt$input)
if (!is.finite(opt$fdr) || opt$fdr <= 0 || opt$fdr > 1) {
  stop("--fdr must be greater than 0 and no greater than 1")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- fread(opt$input)
required_cols <- c("Gene", "pval_column", "ACAT_p")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols)) {
  stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df[pval_column == opt$main_col]
if (!nrow(df)) stop("No rows for the requested analysis")

df[, Gene := as.factor(Gene)]
df[, pval_column := as.factor(pval_column)]

cat("===== CASTIE Step 4 eGene calling start =====\n")

# Main-effect eGene calling
egene_list <- list()

for (col in levels(df$pval_column)) {
  sub <- df[pval_column == col]
  valid <- is.finite(sub$ACAT_p) & sub$ACAT_p > 0 & sub$ACAT_p <= 1
  sub[, qv := NA_real_]
  if (any(valid)) {
    sub[valid, qv := qvalue(ACAT_p, pi0 = 1, lfdr.out = FALSE)$qvalue]
  }
  
  message(col, " qvalue summary:")
  print(summary(sub$qv))
  
  egenes <- sub[!is.na(qv) & qv <= opt$fdr][order(qv)]
  message(col, " eGenes at FDR ", opt$fdr, ": ", nrow(egenes))
  egene_list[[col]] <- as.character(egenes$Gene)
  
  outfile <- file.path(opt$outdir, paste0(opt$prefix, col, "_egene.tsv"))
  fwrite(egenes, outfile, sep = "\t")
  message("Written: ", outfile)
  
  txtfile <- file.path(opt$outdir, paste0(opt$prefix, col, "_egene_genes.txt"))
  writeLines(as.character(egenes$Gene), txtfile)
  message("Written: ", txtfile)
}

message("\nStep 4 complete.")
