# Generic GSEA summary with a custom chromosome-arm analysis
# Purpose: rank genes from differential-expression results and test pathway enrichment
# Tools used: clusterProfiler, ReactomePA, msigdbr, org.Hs.eg.db, bitr

library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)

# Load DE results from a previous edgeR step. Expected columns include gene_symbol,
# logFC, and the unadjusted PValue.
results <- read_tsv("differential_expression_results.tsv")

#### Prepare a signed ranked list ####

# The sign records direction of change; -log10(PValue) weights genes by evidence.
# Use the unadjusted PValue for the ranking metric and BH adjustment in enrichment.
stopifnot(
  all(!is.na(results$gene_symbol)),
  all(is.finite(results$logFC)),
  all(is.finite(results$PValue)),
  all(results$PValue > 0)
)

rank_table <- results %>%
  transmute(
    gene_symbol = gene_symbol,
    rank_score = sign(logFC) * (-log10(PValue))
  ) %>%
  filter(!is.na(gene_symbol), is.finite(rank_score)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# clusterProfiler expects a named numeric vector sorted from highest to lowest.
rank_symbols <- rank_table$rank_score
names(rank_symbols) <- rank_table$gene_symbol
rank_symbols <- sort(rank_symbols, decreasing = TRUE)

#### Standard Hallmark GSEA using gene symbols ####

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks <- GSEA(
  rank_symbols,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
  eps = 0
)

#### Custom chromosome-arm GSEA using MSigDB positional (C1) gene sets ####

# C1 contains positional gene sets. Several source sets can refer to the same
# chromosome arm; normalize their names to a single arm label and pool their
# member genes. For example, source names beginning with CHR1P are combined as
# one "chr1p" term. This creates arm-level unions rather than testing each
# finer-grained positional set independently.
positional_sets <- msigdbr(species = "Homo sapiens", category = "C1") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  filter(!is.na(entrez_gene)) %>%
  mutate(
    source_name_lower = stringr::str_to_lower(gs_name),
    chromosome_arm = stringr::str_replace(
      source_name_lower,
      "^chr([1-9]|1[0-9]|2[0-2])([pq]).*$",
      "chr\\1\\2"
    )
  ) %>%
  # Keep autosomal p- and q-arms only; discard sex chromosomes, mitochondrial
  # sets, and any other C1 sets that do not match the arm naming convention.
  filter(stringr::str_detect(chromosome_arm, "^chr([1-9]|1[0-9]|2[0-2])[pq]$")) %>%
  distinct(chromosome_arm, entrez_gene)

# TERM2GENE maps each custom pathway name to its member Entrez gene IDs.
chromosome_arm_term2gene <- positional_sets %>%
  transmute(term = chromosome_arm, gene = as.character(entrez_gene)) %>%
  distinct(term, gene)

# Convert the symbol-ranked list to Entrez IDs and resolve duplicated mappings.
# Keep one row per symbol and one row per Entrez ID so the named vector has
# unique identifiers, as required by the enrichment procedure.
symbol_to_entrez <- bitr(
  gene = names(rank_symbols),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  distinct(ENTREZID, .keep_all = TRUE)

rank_entrez_table <- symbol_to_entrez %>%
  mutate(rank_score = unname(rank_symbols[SYMBOL])) %>%
  filter(!is.na(rank_score), is.finite(rank_score)) %>%
  arrange(desc(rank_score)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

rank_entrez <- rank_entrez_table$rank_score
names(rank_entrez) <- as.character(rank_entrez_table$ENTREZID)
rank_entrez <- sort(rank_entrez, decreasing = TRUE)

# Larger maxGSSize accommodates broad, merged chromosome-arm sets.
# pvalueCutoff = 1 retains the full result table; interpret adjusted p-values
# when selecting enriched arms. minGSSize = 10 excludes very small sets.
gsea_chromosome_arms <- GSEA(
  rank_entrez,
  TERM2GENE = chromosome_arm_term2gene,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 50000,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
  eps = 0
)

#### Other pathway databases ####

# These use Entrez-ranked genes as well:
# gsea_wikipathways <- gseWP(
#   rank_entrez,
#   organism = "Homo sapiens",
#   pvalueCutoff = 0.2,
#   minGSSize = 10,
#   maxGSSize = 500
# )
# gsea_reactome <- ReactomePA::gsePathway(
#   rank_entrez,
#   pvalueCutoff = 1,
#   pAdjustMethod = "BH",
#   minGSSize = 10,
#   maxGSSize = 500,
#   verbose = TRUE,
#   seed = TRUE,
#   by = "fgsea",
#   eps = 0
# )

# Save the ordinary pathway and custom arm-level enrichment results.
write_tsv(gsea_hallmarks@result, "gsea_hallmarks_results.tsv")
write_tsv(gsea_chromosome_arms@result, "gsea_chromosome_arm_results.tsv")
