############################################################
# MitoCarta × DepMap CERES — essentiality summary per gene
# For each MitoCarta gene: N_lines, Mean/Median/SD CERES,
# p-value (vs 0), FDR; then volcano plot (Mean CERES vs -log10(p))
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# CONFIG: paths (adjust to your files)
# -----------------------------
MITOCARTA_FILE <- "Human.MitoCarta3.0.xlsx"
CERES_FILE     <- "CRISPRGeneEffect.csv"
MODEL_FILE     <- "Model.csv"
OUT_DIR        <- "results"
OUT_TABLE      <- file.path(OUT_DIR, "MitoCarta_CERES_summary.csv")
OUT_VOLCANO    <- file.path(OUT_DIR, "MitoCarta_CERES_volcano.png")

# CERES matrix expected: rows = cell lines (1st col = line ID or name), columns = genes
# If your file is genes × cell lines, set TRANSPOSE_CERES <- TRUE
TRANSPOSE_CERES <- FALSE
USE_CELL_LINE_NAMES <- TRUE   # use Model.csv to set row names to cell line names

dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------
# 1) Load MitoCarta gene list
# -----------------------------
cat("Loading MitoCarta...\n")
mito    <- read_excel(MITOCARTA_FILE, sheet = 1)
mito_genes <- toupper(unique(na.omit(mito$Symbol)))
cat("MitoCarta genes:", length(mito_genes), "\n")

# -----------------------------
# 2) Load CERES matrix
# -----------------------------
cat("Loading CERES matrix...\n")
ceres <- fread(CERES_FILE)
# Assume: 1st column = cell line identifier (DepMap_ID or name), rest = genes
id_col <- names(ceres)[1]
ceres_mat <- as.matrix(ceres[, -1, with = FALSE])
rownames(ceres_mat) <- ceres[[id_col]]
colnames(ceres_mat) <- names(ceres)[-1]

if (TRANSPOSE_CERES) {
  ceres_mat <- t(ceres_mat)
}

# Ensure we have: rows = cell lines, columns = genes
# Heuristic: if more rows than columns, file is genes x lines → transpose
if (nrow(ceres_mat) > ncol(ceres_mat)) {
  ceres_mat <- t(ceres_mat)
  cat("Transposed CERES so that rows = cell lines, columns = genes.\n")
}

# Optional: replace row names with cell line names from Model.csv
if (USE_CELL_LINE_NAMES && file.exists(MODEL_FILE)) {
  model <- fread(MODEL_FILE)
  # Common DepMap Model columns: ModelID, StrippedCellLineName or CellLineName
  id_in_model  <- names(model)[tolower(names(model)) %in% c("modelid", "depmap_id")][1]
  name_in_model <- names(model)[tolower(names(model)) %in% c("strippedcelllinename", "celllinename", "cell_line_name")][1]
  if (is.na(id_in_model)) id_in_model <- names(model)[1]
  if (is.na(name_in_model)) name_in_model <- names(model)[2]
  orig_rows <- rownames(ceres_mat)
  new_names <- model[[name_in_model]][match(orig_rows, model[[id_in_model]])]
  rownames(ceres_mat) <- ifelse(is.na(new_names), orig_rows, new_names)
  cat("Row names set to cell line names from", MODEL_FILE, "\n")
}

# Numeric matrix
storage.mode(ceres_mat) <- "numeric"
cat("CERES dimensions: ", nrow(ceres_mat), "lines x ", ncol(ceres_mat), "genes\n")

# -----------------------------
# 2b) Normalize CERES gene identifiers to HUGO symbols (upper case)
#     This makes them comparable to MitoCarta symbols.
# -----------------------------
genes_raw <- colnames(ceres_mat)

# If columns look like 'GENE (1234)', keep only the gene symbol part
genes_sym <- toupper(gsub("\\s*\\(.*\\)$", "", genes_raw))

# If columns look like Entrez IDs, try mapping Entrez -> HUGO using the same map
only_digits <- grepl("^[0-9]+$", genes_sym)
if (any(only_digits) && file.exists("uniprot_hugo_entrez_id_mapping.csv")) {
  map <- fread("uniprot_hugo_entrez_id_mapping.csv")
  map <- map[!is.na(EntrezID) & !is.na(Symbol)]
  # EntrezID in the mapping is numeric with .0, convert to integer-like character
  map$EntrezID_chr <- as.character(as.integer(map$EntrezID))

  genes_entrez <- genes_sym
  genes_entrez[!only_digits] <- NA_character_

  hugo_from_entrez <- map$Symbol[match(genes_entrez, map$EntrezID_chr)]

  # Where mapping succeeds, replace with HUGO symbols
  replace_idx <- which(!is.na(hugo_from_entrez))
  if (length(replace_idx) > 0) {
    genes_sym[replace_idx] <- toupper(hugo_from_entrez[replace_idx])
  }
}

colnames(ceres_mat) <- genes_sym

# -----------------------------
# 3) Restrict to MitoCarta genes present in CERES
# -----------------------------
genes_in_ceres <- colnames(ceres_mat)
common_genes   <- intersect(mito_genes, genes_in_ceres)
ceres_mito     <- ceres_mat[, common_genes, drop = FALSE]
cat("MitoCarta genes in CERES:", ncol(ceres_mito), "\n")
if (ncol(ceres_mito) == 0) stop("No MitoCarta genes found in CERES. Check gene symbols.")

# -----------------------------
# 4) Per-gene stats: N_lines, Mean, Median, SD, p-value (vs 0), FDR
# -----------------------------
one_sample_p <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  tryCatch(t.test(x, mu = 0)$p.value, error = function(e) NA_real_)
}

stats_list <- lapply(seq_len(ncol(ceres_mito)), function(j) {
  g  <- colnames(ceres_mito)[j]
  v  <- ceres_mito[, j]
  v  <- v[!is.na(v) & is.finite(v)]
  n  <- length(v)
  data.frame(
    Gene         = g,
    N_lines      = n,
    Mean_CERES   = if (n > 0) mean(v) else NA_real_,
    Median_CERES = if (n > 0) median(v) else NA_real_,
    SD_CERES     = if (n > 1) sd(v) else NA_real_,
    p_value      = one_sample_p(ceres_mito[, j]),
    stringsAsFactors = FALSE
  )
})

res <- bind_rows(stats_list)
res$FDR_q <- p.adjust(res$p_value, method = "BH")
res <- res %>% arrange(p_value)

# -----------------------------
# 5) Save summary table
# -----------------------------
fwrite(res, OUT_TABLE)
cat("Saved:", OUT_TABLE, "\n")

# -----------------------------
# 6) Volcano-style plot: X = Mean CERES, Y = -log10(p_value)
# Left = essential genes (negative CERES); Top = significant; Top-left = core dependencies
# -----------------------------
res_plot <- res %>%
  filter(is.finite(p_value) & p_value > 0) %>%
  mutate(
    neglog10p = -log10(p_value),
    label = if_else(FDR_q < 0.05 & Mean_CERES < -0.5, as.character(Gene), "")
  )

p_volcano <- ggplot(res_plot, aes(x = Mean_CERES, y = neglog10p)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "gray40", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "solid") +
  theme_minimal() +
  labs(
    title = "MitoCarta genes: CERES essentiality vs statistical significance",
    subtitle = "Left = essential; Top = significant; Top-left = core mitochondrial dependencies",
    x = "Mean CERES",
    y = "-log10(p-value)"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, box.padding = 0.3)
} else {
  p_volcano <- p_volcano +
    geom_text(aes(label = label), hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE)
}

ggsave(OUT_VOLCANO, plot = p_volcano, width = 9, height = 7, dpi = 300)
cat("Saved:", OUT_VOLCANO, "\n")

cat("Done.\n")
