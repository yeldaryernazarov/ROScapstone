############################################################
# Expression normalization and Essential vs NonEssential score
# Steps:
#   1) log2(TPM + 1)
#   2) Z-score per gene across samples
#   3) Per-sample mean Z for essential_genes and nonessential_genes
#   4) RelativeScore = EssentialScore - NonEssentialScore
# Output: SampleID, EssentialScore, NonEssentialScore, RelativeScore
# Gene symbols = rownames of expression matrix
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# -----------------------------
# CONFIG
# -----------------------------
# Expression matrix file: rows = genes (symbols), columns = samples (SampleID)
# Or: first column = Gene, rest = samples
EXPR_FILE <- "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"

# If expression is already log2(TPM+1), set SKIP_LOG2 = TRUE
SKIP_LOG2 <- TRUE

# Gene lists: paths to files (one gene per line or CSV with 'Gene' column) or character vector
# Default: use MitoCarta as "essential" and a non-essential reference set
ESSENTIAL_GENES_FILE   <- "essential_genes.txt"
NONESSENTIAL_GENES_FILE <- "nonessential_genes.txt"

OUT_DIR  <- "results"
OUT_CSV  <- file.path(OUT_DIR, "Essentiality_RelativeScore.csv")

dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------
# Load gene lists (optional: create from MitoCarta + reference)
# -----------------------------
read_genes <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- fread(path, header = FALSE)
  vec <- as.character(x[[1]])
  vec <- vec[!grepl("^#", vec)]  # drop comment lines
  if (ncol(x) >= 1 && tolower(names(x)[1]) == "gene") return(toupper(unique(na.omit(vec))))
  return(toupper(unique(na.omit(vec))))
}

essential_genes   <- read_genes(ESSENTIAL_GENES_FILE)
nonessential_genes <- read_genes(NONESSENTIAL_GENES_FILE)

# If lists are empty, try MitoCarta for essential and a placeholder for non-essential
if (length(essential_genes) == 0 && file.exists("Human.MitoCarta3.0.xlsx")) {
  suppressPackageStartupMessages(library(readxl))
  mito <- read_excel("Human.MitoCarta3.0.xlsx", sheet = 1)
  essential_genes <- toupper(unique(na.omit(mito$Symbol)))
  cat("Using MitoCarta as essential_genes:", length(essential_genes), "genes\n")
}

if (length(nonessential_genes) == 0) {
  # Common non-essential reference (e.g. Hart non-essential or leave empty → then only EssentialScore)
  nonessential_genes <- character(0)
  cat("No nonessential_genes file found. NonEssentialScore will be 0 for all samples.\n")
}

# -----------------------------
# Load expression matrix (genes = rownames or first column, samples = columns)
# -----------------------------
cat("Loading expression...\n")
expr <- fread(EXPR_FILE)

# Detect layout: if first column name is "Gene" or looks like IDs, use as gene IDs
first_col <- names(expr)[1]
if (tolower(first_col) %in% c("gene", "symbol", "gene_symbol", "geneid")) {
  gene_col <- first_col
  sample_cols <- setdiff(names(expr), gene_col)
  mat <- as.matrix(expr[, ..sample_cols])
  rownames(mat) <- as.character(expr[[gene_col]])
} else {
  # Assume rows are genes (first column is gene ID)
  mat <- as.matrix(expr[, -1, with = FALSE])
  rownames(mat) <- as.character(expr[[1]])
}
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- make.unique(rownames(mat))
cat("Expression: ", nrow(mat), "genes x ", ncol(mat), "samples\n")

# -----------------------------
# 1) log2(TPM + 1) if not already
# -----------------------------
if (!SKIP_LOG2) {
  mat <- log2(mat + 1)
  cat("Applied log2(TPM + 1).\n")
} else {
  cat("Skipping log2 (assuming already transformed).\n")
}

# -----------------------------
# 2) Z-score per gene across samples
# -----------------------------
mat_z <- t(scale(t(mat)))
cat("Z-score normalized per gene across samples.\n")

# -----------------------------
# 3) Per-sample mean Z for essential and nonessential
# -----------------------------
all_genes <- rownames(mat_z)
all_genes_upper <- toupper(all_genes)

essential_idx   <- which(all_genes_upper %in% toupper(essential_genes))
nonessential_idx <- which(all_genes_upper %in% toupper(nonessential_genes))

EssentialScore <- colMeans(mat_z[essential_idx, , drop = FALSE], na.rm = TRUE)
if (length(nonessential_idx) > 0) {
  NonEssentialScore <- colMeans(mat_z[nonessential_idx, , drop = FALSE], na.rm = TRUE)
} else {
  NonEssentialScore <- setNames(rep(0, ncol(mat_z)), colnames(mat_z))
}

# -----------------------------
# 4) RelativeScore = EssentialScore - NonEssentialScore
# -----------------------------
RelativeScore <- EssentialScore - NonEssentialScore

# -----------------------------
# Output dataframe
# -----------------------------
out <- data.frame(
  SampleID         = names(EssentialScore),
  EssentialScore   = EssentialScore,
  NonEssentialScore = NonEssentialScore,
  RelativeScore    = RelativeScore,
  stringsAsFactors = FALSE
)
rownames(out) <- NULL

fwrite(out, OUT_CSV)
cat("Saved:", OUT_CSV, "\n")
cat("Genes used - Essential:", length(essential_idx), " NonEssential:", length(nonessential_idx), "\n")
cat("Done.\n")
