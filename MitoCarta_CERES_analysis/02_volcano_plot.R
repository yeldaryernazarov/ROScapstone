############################################################
# Volcano plot from MitoCarta CERES summary table
# X = Mean CERES, Y = -log10(p-value)
# Run after 01_mitocarta_ceres_summary.R (or use existing summary CSV)
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

SUMMARY_CSV <- "results/MitoCarta_CERES_summary.csv"
OUT_PNG     <- "results/MitoCarta_CERES_volcano.png"

if (!file.exists(SUMMARY_CSV)) {
  stop("Run 01_mitocarta_ceres_summary.R first to create ", SUMMARY_CSV)
}

res <- fread(SUMMARY_CSV)
res_plot <- res %>%
  filter(is.finite(p_value) & p_value > 0) %>%
  mutate(
    neglog10p = -log10(p_value),
    label = if_else(FDR_q < 0.05 & Mean_CERES < -0.5, as.character(Gene), "")
  )

p <- ggplot(res_plot, aes(x = Mean_CERES, y = neglog10p)) +
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
  p <- p + ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, box.padding = 0.3)
} else {
  p <- p + geom_text(aes(label = label), hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE)
}

dir.create("results", showWarnings = FALSE)
ggsave(OUT_PNG, plot = p, width = 9, height = 7, dpi = 300)
cat("Saved:", OUT_PNG, "\n")
