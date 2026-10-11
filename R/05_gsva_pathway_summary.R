# Generic GSVA pathway scoring summary
# Purpose: anonymized example of pathway-level scoring from gene-expression data
# Tools used: GSVA, msigdbr, tidyverse, ggplot2

library(tidyverse)
library(msigdbr)
library(GSVA)

# Read expression matrix: genes x samples, usually TPM or log2-transformed TPM
expr <- read.delim("expression_matrix.tsv", sep = "\t", header = TRUE, row.names = NULL)
expr <- expr %>% column_to_rownames("gene_symbol")

# Remove duplicate genes by keeping the representative row with highest mean expression
expr <- expr %>%
  mutate(mean_expr = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  group_by(gene_symbol) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-mean_expr)

# Build pathway collection from MSigDB Hallmark gene sets
pathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways_list <- split(pathways$gene_symbol, pathways$gs_name)

# Convert to matrix expected by GSVA
expr_matrix <- as.matrix(expr)

# Pathway score calculation
param <- gsvaParam(expr_matrix, pathways_list)
gsva_scores <- gsva(param)

gsva_scores <- t(gsva_scores) %>% as.data.frame()
gsva_scores$sample_id <- rownames(gsva_scores)

# Merge with metadata (purity, group, clinical annotation)
metadata <- read.delim("sample_metadata.tsv", sep = "\t", header = TRUE)
gsva_scores <- merge(gsva_scores, metadata, by = "sample_id", all = FALSE)

# Example correlation test between purity and one pathway score
res <- cor.test(
  gsva_scores$purity,
  gsva_scores$EMT_pathway,
  method = "spearman"
)

# Plot pathway score vs metadata
p <- ggplot(gsva_scores, aes(x = purity, y = EMT_pathway)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "sample purity", y = "GSVA score", title = "EMT pathway") +
  theme_minimal()

# Save outputs
write_tsv(gsva_scores, "gsva_pathway_scores.tsv")
