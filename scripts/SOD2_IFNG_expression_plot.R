############################################################
# SOD2 expression vs interferon-gamma signature
# - Uses existing SOD2 z-scores per patient
# - Combines with interferon-gamma (IFNG) signal, if available
# - Falls back to ESTIMATE ImmuneScore as proxy if IFNG file is absent
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

############################################################
# CONFIG: paths (adjust to your files)
############################################################

# SOD2 per-patient table, already present in the project
# Expected columns (minimum): patient, SOD2 or SOD2_z
SOD2_FILE     <- "результаты/TCGA_Q25Q75_patients_SOD2.tsv"

# Optional: interferon-gamma expression/signature per patient.
# If this file exists, it should contain:
#   - patient
#   - IFNG   (or IFNG_score) — numeric value per patient
#
# Example structure:
# patient,IFNG
# TCGA-XX-XXXX, 1.23
# ...
#
# If this file does NOT exist, the script will instead
# derive an "IFNG proxy" from ESTIMATE ImmuneScore.
IFNG_FILE     <- "результаты/TCGA_IFNG_expression.tsv"

# ESTIMATE scores for all TCGA patients (already in repo)
# Used as a fallback proxy for interferon-gamma activity.
# Expected columns: patient, ImmuneScore (and others)
ESTIMATE_FILE <- "TCGA_ESTIMATE_all_patients.tsv"

# Output directory and figures
OUT_DIR          <- "результаты"
OUT_SCATTER      <- file.path(OUT_DIR, "SOD2_IFNG_scatter.png")
OUT_BOX_IFNG_X   <- file.path(OUT_DIR, "SOD2_by_IFNG_group_boxplot.png")
OUT_BOX_SOD2_X   <- file.path(OUT_DIR, "IFNG_by_SOD2_group_boxplot.png")

dir.create(OUT_DIR, showWarnings = FALSE)

############################################################
# 1) Load SOD2 data
############################################################

if (!file.exists(SOD2_FILE)) {
  stop("SOD2 file not found: ", SOD2_FILE)
}

dt_sod2 <- fread(SOD2_FILE)

if (!"patient" %in% names(dt_sod2)) {
  stop("SOD2 table must contain a 'patient' column.")
}

dt_sod2[, patient := substr(patient, 1, 12)]

# Make sure we have a numeric SOD2_z column
if (!"SOD2_z" %in% names(dt_sod2)) {
  if (!"SOD2" %in% names(dt_sod2)) {
    stop("SOD2 table must contain either 'SOD2_z' or 'SOD2' column.")
  }
  dt_sod2[, SOD2 := as.numeric(SOD2)]
  mu_sod2 <- mean(dt_sod2$SOD2, na.rm = TRUE)
  sd_sod2 <- sd(dt_sod2$SOD2, na.rm = TRUE)
  if (is.na(sd_sod2) || sd_sod2 == 0) {
    stop("Cannot compute z-score for SOD2: zero or NA standard deviation.")
  }
  dt_sod2[, SOD2_z := (SOD2 - mu_sod2) / sd_sod2]
}

dt_sod2 <- dt_sod2[is.finite(SOD2_z)]

############################################################
# 2) Load IFNG / interferon-gamma signal
############################################################

if (file.exists(IFNG_FILE)) {
  # Preferred: direct interferon-gamma signal per patient
  message("Using IFNG expression/signature from: ", IFNG_FILE)
  dt_ifng <- fread(IFNG_FILE)

  if (!"patient" %in% names(dt_ifng)) {
    stop("IFNG table must contain a 'patient' column.")
  }
  dt_ifng[, patient := substr(patient, 1, 12)]

  if ("IFNG" %in% names(dt_ifng)) {
    dt_ifng[, IFNG_value := as.numeric(IFNG)]
  } else if ("IFNG_score" %in% names(dt_ifng)) {
    dt_ifng[, IFNG_value := as.numeric(IFNG_score)]
  } else {
    stop(
      "IFNG table must contain either 'IFNG' or 'IFNG_score' column.\n",
      "Please adjust the code or column names accordingly."
    )
  }

  dt_ifng <- dt_ifng[is.finite(IFNG_value)]

  # Z-score for IFNG
  mu_ifng <- mean(dt_ifng$IFNG_value, na.rm = TRUE)
  sd_ifng <- sd(dt_ifng$IFNG_value, na.rm = TRUE)
  if (is.na(sd_ifng) || sd_ifng == 0) {
    stop("Cannot compute z-score for IFNG: zero or NA standard deviation.")
  }
  dt_ifng[, IFNG_z := (IFNG_value - mu_ifng) / sd_ifng]

  ifng_source_label <- "IFNG expression/signature"

} else {
  # Fallback: use ESTIMATE ImmuneScore as a proxy
  message(
    "IFNG file not found (", IFNG_FILE, ").\n",
    "Falling back to ESTIMATE ImmuneScore as an interferon-gamma proxy."
  )

  if (!file.exists(ESTIMATE_FILE)) {
    stop("ESTIMATE file not found: ", ESTIMATE_FILE)
  }

  est <- fread(ESTIMATE_FILE)
  if (!"patient" %in% names(est)) {
    stop("ESTIMATE table must contain a 'patient' column.")
  }
  est[, patient := substr(patient, 1, 12)]

  if (!"ImmuneScore" %in% names(est)) {
    stop("ESTIMATE table must contain an 'ImmuneScore' column.")
  }

  est[, ImmuneScore := as.numeric(ImmuneScore)]
  est <- est[is.finite(ImmuneScore)]

  mu_ifng <- mean(est$ImmuneScore, na.rm = TRUE)
  sd_ifng <- sd(est$ImmuneScore, na.rm = TRUE)
  if (is.na(sd_ifng) || sd_ifng == 0) {
    stop("Cannot compute z-score for ImmuneScore: zero or NA standard deviation.")
  }

  dt_ifng <- est[, .(patient, IFNG_value = ImmuneScore)]
  dt_ifng[, IFNG_z := (IFNG_value - mu_ifng) / sd_ifng]

  ifng_source_label <- "ImmuneScore (ESTIMATE, IFNG proxy)"
}

dt_ifng <- dt_ifng[is.finite(IFNG_z)]

############################################################
# 3) Merge SOD2 and IFNG
############################################################

dt <- merge(
  dt_sod2[, .(patient, SOD2_z)],
  dt_ifng[, .(patient, IFNG_z)],
  by = "patient",
  all = FALSE
)

dt <- dt[is.finite(SOD2_z) & is.finite(IFNG_z)]

if (nrow(dt) == 0) {
  stop("No overlapping patients with valid SOD2_z and IFNG_z.")
}

message("Merged patients with SOD2 and IFNG signal: ", nrow(dt))

############################################################
# 4) Define IFNG high / low groups (Q25 vs Q75)
############################################################

q25_ifng <- quantile(dt$IFNG_z, 0.25, na.rm = TRUE)
q75_ifng <- quantile(dt$IFNG_z, 0.75, na.rm = TRUE)

dt_groups <- dt[(IFNG_z <= q25_ifng) | (IFNG_z >= q75_ifng)]
dt_groups[, IFNG_group := ifelse(
  IFNG_z >= q75_ifng,
  "IFNG_high_25pct",
  "IFNG_low_25pct"
)]

message("IFNG group sizes:")
print(table(dt_groups$IFNG_group))

############################################################
# 5) Scatter plot: SOD2_z vs IFNG_z
############################################################

p_scatter <- ggplot(dt, aes(x = IFNG_z, y = SOD2_z)) +
  geom_point(alpha = 0.6, size = 1.8, color = "#2C7BB6") +
  geom_smooth(method = "lm", se = TRUE, color = "#D7191C") +
  theme_minimal() +
  labs(
    title    = "SOD2 expression vs interferon-gamma signal",
    subtitle = ifng_source_label,
    x        = "IFNG z-score",
    y        = "SOD2 z-score"
  )

ggsave(OUT_SCATTER, plot = p_scatter, width = 7, height = 6, dpi = 300)
cat("Saved:", OUT_SCATTER, "\n")

############################################################
# 6) Boxplot 1: SOD2_z by IFNG high/low (Q25 vs Q75)
############################################################

p_box_ifng_x <- ggplot(dt_groups, aes(x = IFNG_group, y = SOD2_z, fill = IFNG_group)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  scale_fill_manual(values = c("IFNG_low_25pct" = "#2C7BB6",
                               "IFNG_high_25pct" = "#D7191C")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title    = "SOD2 expression by interferon-gamma high vs low (Q25 vs Q75)",
    subtitle = ifng_source_label,
    x        = "",
    y        = "SOD2 z-score"
  )

ggsave(OUT_BOX_IFNG_X, plot = p_box_ifng_x, width = 6, height = 5, dpi = 300)
cat("Saved:", OUT_BOX_IFNG_X, "\n")

############################################################
# 7) Boxplot 2: IFNG_z by SOD2 high/low (Q25 vs Q75)
############################################################

# Define SOD2 high / low (Q25 vs Q75)
q25_sod2 <- quantile(dt$SOD2_z, 0.25, na.rm = TRUE)
q75_sod2 <- quantile(dt$SOD2_z, 0.75, na.rm = TRUE)

dt_sod2_groups <- dt[(SOD2_z <= q25_sod2) | (SOD2_z >= q75_sod2)]
dt_sod2_groups[, SOD2_group := ifelse(
  SOD2_z >= q75_sod2,
  "SOD2_high_25pct",
  "SOD2_low_25pct"
)]

message("SOD2 group sizes:")
print(table(dt_sod2_groups$SOD2_group))

p_box_sod2_x <- ggplot(dt_sod2_groups, aes(x = SOD2_group, y = IFNG_z, fill = SOD2_group)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  scale_fill_manual(values = c("SOD2_low_25pct" = "#2C7BB6",
                               "SOD2_high_25pct" = "#D7191C")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title    = "Interferon-gamma signal by SOD2 high vs low (Q25 vs Q75)",
    subtitle = ifng_source_label,
    x        = "",
    y        = "IFNG z-score"
  )

ggsave(OUT_BOX_SOD2_X, plot = p_box_sod2_x, width = 6, height = 5, dpi = 300)
cat("Saved:", OUT_BOX_SOD2_X, "\n")

############################################################
# DONE
############################################################
cat("Done.\n")

