# PSG expression analysis in Preeclampsia (Control vs PE)

This repository contains lightweight, reproducible scripts to quantify and visualize **Pregnancy-Specific Glycoprotein (PSG)** gene expression in **preeclampsia (PE)** versus **control** placenta RNA-seq expression tables (TPM). The main output is a set of per-gene plots and summary statistics using a non-parametric test (Wilcoxon rank-sum), plus exported matrices for downstream figure-making.

## What this does

Given:
- a TPM matrix (`Gene TPM.csv`) and
- an Ensembl ID ↔ gene symbol mapping table (`ensembl_human_ID.txt`)

the pipeline:
1. maps Ensembl IDs to gene symbols (when applicable),
2. extracts all genes matching `^PSG`,
3. transforms values to `log2(TPM+1)`,
4. generates:
   - one multi-panel PDF containing all PSG genes,
   - individual PNG plots per PSG gene,
5. exports:
   - raw PSG TPM matrix,
   - log2(TPM+1) PSG matrix,
   - row-centered PSG matrix (per-gene mean subtracted),
   - stats table (means, delta, Wilcoxon p-value + BH adjusted p-value).

## Repository structure

├── README.md
├── scripts/
│ └── psg_pe_expression.R
└── PSG_PE_plots/ # created automatically after running
├── PSG_all_boxplots.pdf
├── PSG*_Expression_Control_PE.png
├── PSG_TPM_raw.tsv
├── PSG_log2TPM1.tsv
├── PSG_log2TPM1_rowCentered.tsv
└── PSG_stats_Control_vs_PE.tsv


## Input requirements

### 1) TPM table (`Gene TPM.csv`)
A CSV with:
- first column = gene identifier (Ensembl ID **or** gene symbol)
- remaining columns = samples (numeric TPM values)

Example (conceptual):

### 2) Ensembl mapping (`ensembl_human_ID.txt`)
A tab-delimited file with at least:
- one column containing Ensembl IDs (e.g., `ENSG...`)
- one column containing gene symbols (e.g., `PSG9`)

The script auto-detects common column names (`ensembl`, `gene_id`, `symbol`, `gene_name`, etc.).  
If your mapping file has a different format, adjust the `ens_id_col` / `sym_col` selection in the script.

## Installation

### R version
Tested on R ≥ 4.0 (should work on newer versions).

### R packages
Only standard packages are used:
- `dplyr`
- `stringr`

Install if needed:
```r
install.packages(c("dplyr", "stringr"))

Run

Put your input files somewhere accessible (or keep them in your Dropbox path).

Edit the top of the script to point to the correct paths:

tpm_file <- "path/to/Gene TPM.csv"
ens_file <- "path/to/ensembl_human_ID.txt"


Run:

Rscript scripts/psg_pe_expression.R


Outputs will appear in:

PSG_PE_plots/

Group definition (Control vs PE)

By default, the script assumes:

first 8 columns are Control

remaining columns are PE

If your sample order is different, edit:

n_ctrl <- 8
ctrl_cols <- 1:n_ctrl
pe_cols <- (n_ctrl + 1):ncol(psg_log2)


(or define ctrl_cols/pe_cols by sample names).

Statistics

Per-gene test: Wilcoxon rank-sum (robust for non-normal TPM distributions).

Multiple testing correction: BH (FDR) across PSG genes.

The results table is written to:
PSG_PE_plots/PSG_stats_Control_vs_PE.tsv

Notes / caveats

TPM values are not counts, so this script focuses on expression comparison and visualization, not DESeq2-style modeling.

If you have raw counts, use a count-based framework (DESeq2/edgeR) for differential expression.

The script silently converts non-numeric entries to 0 after coercion—if your table contains non-expression columns, remove them first or keep only numeric sample columns.

Citation

If you use this code in a manuscript or preprint, please cite the dataset/source used for the TPM matrix and include a short methods note such as:

PSG genes were extracted from TPM tables, log2(TPM+1) transformed, and compared between Control and PE groups using Wilcoxon rank-sum with BH correction.

Contact

Maintained by: Manvendra Singh
