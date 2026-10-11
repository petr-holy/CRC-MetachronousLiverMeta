#!/usr/bin/env bash
# Reporting summary only; replace generic placeholders with study-specific inputs.

THREADS=number_of_threads
REFERENCE=reference_GRCh38.fa
TARGET_REGIONS=target_regions.bed
ACCESS_REGIONS=accessible_regions.bed
GENE_ANNOTATION=refFlat.txt
ANTITARGET_MASK=excluded_regions.bed
COHORT_BAMS=(sample_1.bam sample_2.bam)
REFERENCE_BAMS=(reference_sample_1.bam reference_sample_2.bam)

# Prepare target intervals from the design BED file.
tail -n +3 "$TARGET_REGIONS" > cnv_targets.bed

# Estimate copy number using CNVkit with a flat reference.
cnvkit.py batch \
  --processes "$THREADS" \
  --fasta "$REFERENCE" \
  --targets cnv_targets.bed \
  --annotate "$GENE_ANNOTATION" \
  --short-names \
  --antitargets "$ANTITARGET_MASK" \
  --output-dir cnv_flat_reference/ \
  --output-reference flat_reference.cnn \
  --drop-low-coverage \
  --segment-method cbs \
  --access "$ACCESS_REGIONS" \
  --scatter --diagram \
  "${COHORT_BAMS[@]}" --normal

# Review coverage and segmentation metrics for the flat-reference run.
cnvkit.py metrics cnv_flat_reference/*.cnr \
  -s cnv_flat_reference/*.cns

# Build a pooled reference from selected samples.
cnvkit.py batch \
  --processes "$THREADS" \
  --fasta "$REFERENCE" \
  --targets cnv_targets.bed \
  --annotate "$GENE_ANNOTATION" \
  --short-names \
  --output-dir cnv_pooled_reference/ \
  --output-reference pooled_reference.cnn \
  --drop-low-coverage \
  --segment-method cbs \
  --access "$ACCESS_REGIONS" \
  --normal "${REFERENCE_BAMS[@]}"

# Call copy-number segments in cohort samples against the pooled reference.
cnvkit.py batch \
  --processes "$THREADS" \
  --output-dir cnv_pooled_calls/ \
  --reference pooled_reference.cnn \
  --drop-low-coverage \
  --segment-method cbs \
  --scatter --diagram \
  "${COHORT_BAMS[@]}"
