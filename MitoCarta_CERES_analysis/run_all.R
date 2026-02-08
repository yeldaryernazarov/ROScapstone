############################################################
# Run MitoCarta CERES summary + volcano (scripts 01 and 02)
# Usage: from MitoCarta_CERES_analysis dir, run
#   Rscript run_all.R
############################################################

script_dir <- if (exists("script_dir")) script_dir else getwd()
setwd(script_dir)

cat("Running 01_mitocarta_ceres_summary.R...\n")
source("01_mitocarta_ceres_summary.R", local = new.env())

cat("Running 02_volcano_plot.R...\n")
source("02_volcano_plot.R", local = new.env())

cat("All done. Check results/ folder.\n")
