#!/usr/bin/env bash

# Set generic inputs and output locations.
SAMPLE_ID="<sample_id>"
FASTQ_R1="<input_R1.fastq.gz>"
FASTQ_R2="<input_R2.fastq.gz>"
TRIM_DIR="<trimmed_reads_dir>"
WORK_DIR="<alignment_output_dir>"
REFERENCE="<reference_genome.fasta>"
ADAPTERS="<paired_end_adapter_fasta>"
TARGET_BED="<target_regions.bed>"
BAIT_BED="<bait_regions.bed>"
INTERVALS="<analysis_intervals.bed>"
KNOWN_SITE_1="<known_sites_1.vcf.gz>"
KNOWN_SITE_2="<known_sites_2.vcf.gz>"
KNOWN_SITE_3="<known_sites_3.vcf.gz>"
KNOWN_SITE_4="<known_sites_4.vcf.gz>"
KNOWN_SITE_5="<known_sites_5.vcf.gz>"


# Convert capture design BED files for hybrid-selection metrics.
gatk BedToIntervalList --INPUT "$TARGET_BED" \
  --SEQUENCE_DICTIONARY "<reference_sequence_dictionary.dict>" \
  --OUTPUT "<target_intervals.txt>"
gatk BedToIntervalList --INPUT "$BAIT_BED" \
  --SEQUENCE_DICTIONARY "<reference_sequence_dictionary.dict>" \
  --OUTPUT "<bait_intervals.txt>"

# Collect alignment and hybrid-selection metrics.
gatk CollectMultipleMetrics --VALIDATION_STRINGENCY LENIENT \
  --INPUT "${WORK_DIR}/${SAMPLE_ID}_analysis_ready.bam" \
  --OUTPUT "${WORK_DIR}/${SAMPLE_ID}_metrics" \
  --REFERENCE_SEQUENCE "$REFERENCE" --FILE_EXTENSION .txt
gatk CollectHsMetrics --BAIT_INTERVALS "<bait_intervals.txt>" \
  --TARGET_INTERVALS "<target_intervals.txt>" \
  --INPUT "${WORK_DIR}/${SAMPLE_ID}_analysis_ready.bam" \
  --OUTPUT "${WORK_DIR}/${SAMPLE_ID}_hybrid_selection_metrics.txt" \
  --METRIC_ACCUMULATION_LEVEL ALL_READS \
  --REFERENCE_SEQUENCE "$REFERENCE" --VALIDATION_STRINGENCY LENIENT

# Summarize target-region depth at selected coverage thresholds.
gatk DepthOfCoverage --intervals "$TARGET_BED" \
  -I "${WORK_DIR}/${SAMPLE_ID}_analysis_ready.bam" \
  -O "${WORK_DIR}/${SAMPLE_ID}_depth_of_coverage.txt" \
  --reference "$REFERENCE" --seconds-between-progress-updates 10 \
  --omit-interval-statistics true --omit-locus-table true \
  --omit-depth-output-at-each-base true \
  --summary-coverage-threshold 15 --summary-coverage-threshold 30 \
  --summary-coverage-threshold 50 --summary-coverage-threshold 75 \
  --summary-coverage-threshold 100 --summary-coverage-threshold 200

# Combine per-sample coverage summaries.
head -n 1 "<representative_sample>.sample_summary" > "<combined_coverage.csv>"
head -n 2 -q *sample_summary | grep -v "above" >> "<combined_coverage.csv>"
