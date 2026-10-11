# Generic RNA-seq differential expression summary
# Purpose: anonymized example of edgeR preprocessing, filtering, and fitting multiple variables
# Tools used: edgeR, EnhancedVolcano, tidyverse

library(edgeR)
library(EnhancedVolcano)
library(tidyverse)

# Read count matrix; rows = genes, columns = samples
count_matrix <- read.delim("counts_matrix.tsv", sep = "\t", header = TRUE, row.names = NULL)
count_matrix <- count_matrix %>% column_to_rownames("gene_id")

# Sample metadata table: one row per sample, includes group labels and covariates
metadata <- read.delim("sample_metadata.tsv", sep = "\t", header = TRUE)

# Example grouping variable
sample_group <- factor(metadata$group_label)

# Create DGEList object
obj <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = sample_group)

# Preprocessing: filter low-expression genes
keep <- filterByExpr(obj)
obj <- obj[keep, , keep.lib.sizes = FALSE]

# Normalization and dispersion estimation
obj <- calcNormFactors(obj, method = "TMM")
obj <- estimateDisp(obj)

# Simple two-group comparison
result_2group <- exactTest(obj, pair = c("group_0", "group_1"))
result_table <- topTags(
  result_2group,
  n = nrow(result_2group$table),
  adjust.method = "BH",
  sort.by = "PValue"
)$table

# Multi-variable modelling example
# This mirrors scripts that loop over several covariates and fit a model per variable
variables <- c("variable_a", "variable_b", "variable_c")
for (var in variables) {
  design <- model.matrix(~ metadata[[var]] + additional_covariate)
  fit <- glmQLFit(obj, design)
  contrast <- makeContrasts(variable_a_vs_reference = metadata[[var]] - 0, levels = design)
  qlf <- glmQLFTest(fit, contrast = contrast)
  de_table <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
}

# Volcano plot export
volcano_plot <- EnhancedVolcano(
  de_table,
  lab = rownames(de_table),
  x = "logFC",
  y = "FDR",
  title = "group comparison",
  subtitle = NULL,
  ylab = "-log10 adjusted p-value",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  arrowheads = FALSE,
  legendLabels = c("NS", "Log2FC>1", "Significant", "Significant & Log2FC>1"),
  legendPosition = "bottom",
  colAlpha = 0.9,
  drawConnectors = TRUE,
  boxedLabels = TRUE
)

# Export results
write_tsv(result_table, "differential_expression_results.tsv")
