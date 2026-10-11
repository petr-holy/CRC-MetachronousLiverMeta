# Generic RNA-seq deconvolution summary
# Purpose: anonymized example of estimating immune/stromal content from bulk gene expression
# Tools used: immunedeconv, EPIC, tidyverse

library(tidyverse)
library(immunedeconv)

# Read TPM-like expression matrix for the bulk samples
expr <- read_tsv("expression_matrix.tsv")
expr <- expr %>% column_to_rownames("gene_symbol")

# Select the sample subset to analyze (e.g., tumor-enriched or high-purity samples)
selected_samples <- colnames(expr)[grepl("sample_group|high_quality", colnames(expr))]
expr_subset <- expr[, selected_samples, drop = FALSE]

# Optional: remove duplicate gene IDs by choosing the representative value for each gene
expr_subset <- expr_subset[!duplicated(rownames(expr_subset)), , drop = FALSE]

# Deconvolution via immunedeconv wrapper
# EPIC: method = "epic", tumor = TRUE
# MCP-counter: method = "mcp_counter", log_transform = TRUE
# ESTIMATE: method = "estimate", log_transform = TRUE

deconv_epic <- deconvolute(
  expr_subset,
  method = "epic",
  tumor = TRUE
)

deconv_mcpcounter <- deconvolute(
  expr_subset,
  method = "mcp_counter",
  log_transform = TRUE,
  tumor = TRUE
)

deconv_estimate <- deconvolute(
  expr_subset,
  method = "estimate",
  log_transform = TRUE,
  tumor = TRUE
)

# Alternative direct EPIC call
# library(EPIC)
# epic_fit <- EPIC(expr_subset)

# Save outputs
write_tsv(deconv_epic, "deconvolution_epic.tsv")
write_tsv(deconv_mcpcounter, "deconvolution_mcp_counter.tsv")
write_tsv(deconv_estimate, "deconvolution_estimate.tsv")

# Typical downstream workflow:
# - join deconvolution estimates with purity or pathway scores
# - compare estimated immune/stromal fractions against GSVA pathway-level scores
# - visualize correlations with ggplot2 or scatterplots
