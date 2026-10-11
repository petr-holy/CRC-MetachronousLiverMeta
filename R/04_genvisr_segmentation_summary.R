# Generic GenVisR segmentation visualization summary
# Purpose: anonymized examples of copy-number frequency and segment-level plots
# Tools used: GenVisR, tidyverse, scales, svglite

library(tidyverse)
library(GenVisR)
library(scales)
library(svglite)

# SEG-like input typically includes sample, chromosome, start, end, and segmean.
seg_data <- read.table("segmented_copynumber.tsv", header = TRUE, sep = "\t")

# Optionally subset samples before plotting (e.g., by cohort or quality group).
# Replace the predicate with a real, non-sensitive selection rule as needed.
selected_samples <- unique(seg_data$sample)
seg_filtered <- seg_data %>% dplyr::filter(sample %in% selected_samples)

chromosomes <- paste0("chr", 1:22)

#### cnFreq: frequency of gains and losses across samples ####

# Thresholds around zero define the neutral copy-number interval.
# Repeat at multiple cutoffs to compare how the definition affects frequencies.
for (cutoff in c(0.1, 0.2, 0.3)) {
  png(
    paste0("cn_freq_cutoff_", cutoff, ".png"),
    width = 42, height = 24, units = "cm", res = 300
  )
  GenVisR::cnFreq(
    seg_filtered,
    genome = "hg38",
    CN_low_cutoff = -cutoff,
    CN_high_cutoff = cutoff,
    plotChr = chromosomes
  )
  dev.off()
}

#### cnSpec: segment-level copy-number profile by sample ####

# CNscale = "relative" displays segment values relative to the neutral baseline.
# Use genome = "hg38" when the input chromosome naming/coordinates are compatible.
svglite::svglite("cn_spec_all_chromosomes.svg", width = 16, height = 12)
GenVisR::cnSpec(
  seg_filtered,
  genome = "hg38",
  CNscale = "relative"
)
dev.off()

# Restrict to selected chromosomes for a more focused profile plot.
# The example data uses chromosome labels without the "chr" prefix in the
# segmentation table; adapt this filter if your input uses a different format.
selected_chromosomes <- c("3", "5", "16")
seg_selected_chromosomes <- seg_filtered %>%
  dplyr::filter(chromosome %in% selected_chromosomes)

# Optional custom genome boundary table (one row per plotted chromosome):
# chromosome contains labels such as "chr3"; start/end are coordinate bounds.
# Supplying y allows the plot to use explicit chromosome boundaries.
genome_boundaries <- read.delim("genome_boundaries.tsv", sep = "\t")
genome_boundaries_selected <- genome_boundaries %>%
  dplyr::filter(chromosome %in% paste0("chr", selected_chromosomes))

svglite::svglite("cn_spec_selected_chromosomes.svg", width = 8, height = 12)
GenVisR::cnSpec(
  seg_selected_chromosomes,
  y = genome_boundaries_selected,
  CNscale = "relative",
  plotLayer = scale_fill_gradientn(
    "Log2 ratio",
    colours = c("blue4", "white", "red4"),
    values = scales::rescale(c(-2.8, 0, 2.8)),
    limits = c(-0.8, 0.8),
    oob = scales::squish
  )
)
dev.off()

# The plotLayer example uses a diverging color scale and squishes values outside
# the specified limits to the endpoint colors. Omit plotLayer to use defaults.
