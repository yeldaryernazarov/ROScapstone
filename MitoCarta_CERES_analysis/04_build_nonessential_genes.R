############################################################
# Build nonessential_genes.txt from MitoCarta_CERES_summary
# - Reads results/MitoCarta_CERES_summary.csv
# - Selects genes with Median_CERES >= 0
# - Writes one gene symbol per line to nonessential_genes.txt
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# CONFIG
# -----------------------------
SUMMARY_FILE <- file.path("results", "MitoCarta_CERES_summary.csv")
OUT_FILE     <- "nonessential_genes.txt"

# -----------------------------
# 1) Load summary table
# -----------------------------
if (!file.exists(SUMMARY_FILE)) {
  stop("Summary file not found: ", SUMMARY_FILE,
       "\nRun 01_mitocarta_ceres_summary.R first, or update SUMMARY_FILE.")
}

cat("Loading summary table from:", SUMMARY_FILE, "\n")
dt <- fread(SUMMARY_FILE)

required_cols <- c("Gene", "Median_CERES")
missing_cols  <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in summary table: ",
       paste(missing_cols, collapse = ", "))
}

# -----------------------------
# 2) Select non-essential candidates (Median_CERES >= 0)
# -----------------------------
dt[, Gene := toupper(Gene)]

nonessential_genes <- dt[
  !is.na(Median_CERES) & is.finite(Median_CERES) & Median_CERES >= 0,
  unique(Gene)
]

cat("Genes with Median_CERES >= 0:", length(nonessential_genes), "\n")

# -----------------------------
# 3) Write nonessential_genes.txt
# -----------------------------
con <- file(OUT_FILE, open = "w", encoding = "UTF-8")
on.exit(close(con), add = TRUE)

header_lines <- c(
  "# One gene symbol per line. Reference non-essential genes (e.g. from Hart et al. or DepMap).",
  "# Auto-generated from results/MitoCarta_CERES_summary.csv using Median_CERES >= 0.",
  "#"
)

writeLines(header_lines, con)

if (length(nonessential_genes) > 0) {
  writeLines(nonessential_genes, con)
}

cat("Wrote", length(nonessential_genes), "genes to", OUT_FILE, "\n")

cat("Done.\n")

