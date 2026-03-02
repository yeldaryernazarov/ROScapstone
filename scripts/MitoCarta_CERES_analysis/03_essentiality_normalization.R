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
# Load expression matrix and map genes to HUGO symbols (like ssGSEA script)
# -----------------------------
cat("Loading expression...\n")
expr <- fread(EXPR_FILE)
cat("Expression dimensions (rows x cols):", dim(expr), "\n")

# In OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv:
#  - Rows are samples
#  - First 4 columns are metadata
#  - Remaining columns are genes like 'GENE (1234)' where 1234 is EntrezID
expr_mat <- as.matrix(expr[, -c(1:4), with = FALSE])
expr_mat <- apply(expr_mat, 2, function(x) as.numeric(as.character(x)))
cat("Converted expr to numeric matrix. NA count:", sum(is.na(expr_mat)), "\n")

# Use SampleID (second column) as sample identifier if present, otherwise first column
sample_id_col <- if ("SampleID" %in% names(expr)) {
  "SampleID"
} else {
  names(expr)[2]
}
rownames(expr_mat) <- as.character(expr[[sample_id_col]])

# Parse Entrez IDs from column names and map to HUGO, exactly as in ssGSEA pipeline
gene_names <- colnames(expr_mat)
entrez_ids <- sub(".*\\((\\d+)\\)$", "\\1", gene_names)

if (file.exists("uniprot_hugo_entrez_id_mapping.csv")) {
  map <- fread("uniprot_hugo_entrez_id_mapping.csv")
  map <- map[!is.na(EntrezID) & !is.na(Symbol)]
  hugo_symbols <- map$Symbol[match(entrez_ids, map$EntrezID)]
  hugo_symbols <- toupper(hugo_symbols)
  colnames(expr_mat) <- hugo_symbols
} else {
  # Fallback: strip trailing '(1234)' and upper-case
  colnames(expr_mat) <- toupper(gsub("\\s*\\(.*\\)$", "", gene_names))
}

# Drop unmapped / empty gene columns
valid_cols <- !is.na(colnames(expr_mat)) & colnames(expr_mat) != ""
expr_mat <- expr_mat[, valid_cols, drop = FALSE]

# Drop all-zero genes (never expressed)
nonzero_cols <- colSums(expr_mat, na.rm = TRUE) > 0
expr_mat <- expr_mat[, nonzero_cols, drop = FALSE]
cat("Columns after mapping and filtering:", ncol(expr_mat), "\n")

# Transpose to genes x samples for downstream per-gene Z-scores
mat <- t(expr_mat)
rownames(mat) <- make.unique(rownames(mat))
cat("Expression matrix (genes x samples):", nrow(mat), "x", ncol(mat), "\n")

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
