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
THREADS="<number_of_threads>"
KNOWN_SITE_1="<known_sites_1.vcf.gz>"
KNOWN_SITE_2="<known_sites_2.vcf.gz>"
KNOWN_SITE_3="<known_sites_3.vcf.gz>"
KNOWN_SITE_4="<known_sites_4.vcf.gz>"
KNOWN_SITE_5="<known_sites_5.vcf.gz>"

# Trim paired-end reads with AGeNT.
AGENT="<agent_executable>"
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_ID}_trimmed_R1.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_ID}_trimmed_R2.fastq.gz"

mkdir -p "$TRIM_DIR"
"$AGENT" trim -v2 -polyG 10 \
  -fq1 "$FASTQ_R1" -fq2 "$FASTQ_R2" \
  -out "${TRIM_DIR}/${SAMPLE_ID}_trimmed"

# Run read-level quality control.
mkdir -p "${TRIM_DIR}/fastqc"
fastqc -t "$THREADS" -o "${TRIM_DIR}/fastqc" "$TRIMMED_R1" "$TRIMMED_R2"

# Align reads, repair mate information, sort, and index.
bwa mem -t "$THREADS" -M \
  -R "@RG\tID:<read_group_id>\tPL:ILLUMINA\tLB:<library_id>\tSM:${SAMPLE_ID}" \
  "$REFERENCE" \
  "$TRIMMED_R1" \
  "$TRIMMED_R2" |
  samtools view -b --threads "$THREADS" - |
  samtools fixmate -m --threads "$THREADS" --output-fmt BAM - - |
  samtools sort --threads "$THREADS" --output-fmt BAM \
    -T "${WORK_DIR}/${SAMPLE_ID}_sort_tmp" \
    -o "${WORK_DIR}/${SAMPLE_ID}_sorted.bam"
samtools index -@ "$THREADS" "${WORK_DIR}/${SAMPLE_ID}_sorted.bam"

# Deduplicate and generate consensus reads with AGeNT CReaK.
"$AGENT" creak \
  --input-read-filtering \
  --consensus-read-filtering \
  --consensus-mode HYBRID \
  --MBC-mismatch 1 \
  --bed-file "$INTERVALS" \
  --output-bam-file "${WORK_DIR}/${SAMPLE_ID}_creak.bam" \
  "${WORK_DIR}/${SAMPLE_ID}_sorted.bam"

# Calculate and apply base-quality score recalibration.
gatk BaseRecalibrator --intervals "$INTERVALS" --reference "$REFERENCE" \
  --known-sites "$KNOWN_SITE_1" --known-sites "$KNOWN_SITE_2" \
  --known-sites "$KNOWN_SITE_3" --known-sites "$KNOWN_SITE_4" \
  --known-sites "$KNOWN_SITE_5" \
  --input "${WORK_DIR}/${SAMPLE_ID}_creak.bam" \
  --output "${WORK_DIR}/${SAMPLE_ID}_recalibration_table"
gatk ApplyBQSR --intervals "$INTERVALS" --reference "$REFERENCE" \
  --input "${WORK_DIR}/${SAMPLE_ID}_creak.bam" \
  --bqsr-recal-file "${WORK_DIR}/${SAMPLE_ID}_recalibration_table" \
  --output "${WORK_DIR}/${SAMPLE_ID}_analysis_ready.bam"
