# Generic RNA-seq PCA and filtering summary
# Purpose: anonymized example of expression filtering, log transformation, and PCA clustering
# Tools used: tidyverse, base R, prcomp

library(tidyverse)

# Read counts matrix; columns are samples, rows are genes
count_matrix <- read.delim("counts_matrix.tsv", sep = "\t", header = TRUE, row.names = NULL)
count_matrix <- count_matrix %>% column_to_rownames("gene_id")

# Filter out low-quality or excluded samples
bad_samples <- c("sample_01", "sample_02")
count_matrix <- count_matrix[, !colnames(count_matrix) %in% bad_samples]

# Keep genes with at least 10 reads in at least 10% of samples
min_reads <- 10
min_fraction <- 0.10
keep_genes <- rowSums(count_matrix >= min_reads) >= (ncol(count_matrix) * min_fraction)
count_matrix <- count_matrix[keep_genes, ]

# Optional: keep only a subset of sample groups for clustering
expr_matrix <- count_matrix[, grepl("group_a|group_b", colnames(count_matrix))]

# Transform for PCA comparability
expr_log <- log2(expr_matrix + 0.001)

# PCA on transposed matrix: rows = samples, columns = genes
samples_matrix <- t(expr_log)
pca_fit <- prcomp(samples_matrix, center = TRUE, scale. = TRUE)

summary(pca_fit)
plot(pca_fit$x[, 1], pca_fit$x[, 2], pch = 16, col = "steelblue")

# Alternative quick check: no scaling if desired
pca_fit_nonscaled <- prcomp(samples_matrix, center = TRUE, scale. = FALSE)

# Typical plotting workflow from RNA-seq clustering scripts:
# - remove outlier samples
# - inspect PC1/PC2 and PC3/PC4
# - color by grouping metadata
# - export PNG/SVG outputs for review
