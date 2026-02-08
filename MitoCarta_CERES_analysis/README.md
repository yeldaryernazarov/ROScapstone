# MitoCarta × DepMap CERES analysis

Scripts to compute **essentiality of MitoCarta genes** across DepMap cell lines using CERES scores, and to run **expression-based essentiality normalization**.

## Data requirements

Place in this folder (or set paths in each script):

| File | Description |
|------|-------------|
| `Human.MitoCarta3.0.xlsx` | MitoCarta 3.0 gene list (sheet 1, column `Symbol`) |
| `CRISPRGeneEffect.csv` | DepMap CERES matrix: **rows = cell lines** (1st col = DepMap ID), **columns = genes** |
| `Model.csv` | DepMap model metadata (columns e.g. `ModelID`, `StrippedCellLineName`) to map ID → cell line name |
| (Optional) `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | Expression matrix for normalization script (genes × samples) |

DepMap data: [DepMap Portal](https://depmap.org/portal/download/) (e.g. CRISPR Gene Effect, Model).

## Scripts

### 1. `01_mitocarta_ceres_summary.R`

- Loads MitoCarta genes (~1500) and DepMap CERES matrix (~2000 lines × genes).
- Optionally maps row names to **cell line names** using `Model.csv`.
- For each MitoCarta gene present in CERES:
  - **N_lines** — number of cell lines with non-NA CERES
  - **Mean_CERES**, **Median_CERES**, **SD_CERES**
  - **p_value** — one-sample test (CERES ≠ 0 across lines)
  - **FDR_q** — Benjamini–Hochberg across MitoCarta genes
- Writes **`results/MitoCarta_CERES_summary.csv`**.
- Builds **volcano plot**: X = Mean CERES, Y = −log10(p-value)  
  - Left = essential genes (negative CERES)  
  - Top = statistically significant  
  - Top-left = core mitochondrial dependencies  
- Saves **`results/MitoCarta_CERES_volcano.png`**.

**Config (top of script):** `MITOCARTA_FILE`, `CERES_FILE`, `MODEL_FILE`, `TRANSPOSE_CERES`, `USE_CELL_LINE_NAMES`.

### 2. `02_volcano_plot.R`

Re-builds the volcano plot from an existing `results/MitoCarta_CERES_summary.csv` (e.g. after editing or filtering). Run after `01_mitocarta_ceres_summary.R` or point `SUMMARY_CSV` to your table.

### 3. `03_essentiality_normalization.R`

Expression-based scoring:

1. **log2(TPM + 1)** (optional if already transformed; set `SKIP_LOG2 = TRUE`).
2. **Z-score** each gene across samples.
3. Per-sample **mean Z-score** for:
   - `essential_genes`
   - `nonessential_genes`
4. **RelativeScore** = EssentialScore − NonEssentialScore.

**Output:** `results/Essentiality_RelativeScore.csv` with columns:  
`SampleID`, `EssentialScore`, `NonEssentialScore`, `RelativeScore`.

**Gene lists:**  
- `essential_genes.txt` — one gene symbol per line (or CSV with column `Gene`).  
- `nonessential_genes.txt` — same format.  
If `essential_genes.txt` is missing and `Human.MitoCarta3.0.xlsx` exists, MitoCarta is used as the essential set.  
If `nonessential_genes.txt` is missing, NonEssentialScore is set to 0.

**Expression matrix:** Genes as row names (or first column `Gene`), samples as columns. Set `EXPR_FILE` and `SKIP_LOG2` at the top of the script.

## Run order

```bash
cd MitoCarta_CERES_analysis
Rscript 01_mitocarta_ceres_summary.R
Rscript 02_volcano_plot.R              # optional
Rscript 03_essentiality_normalization.R # optional; needs expr matrix + gene lists
```

## R packages

- **01, 02:** `data.table`, `readxl`, `dplyr`, `ggplot2`; optional `ggrepel` for volcano labels.
- **03:** `data.table`, `dplyr`, `tidyr`; `readxl` only if using MitoCarta as essential list.

Install optional label repulsion:

```r
install.packages("ggrepel")
```
