#!/usr/bin/env bash

# Demonstration only; replace placeholders with real paths before use.
# Example assumes paired-end FASTQ files.

SAMPLE_ID="<sample_id>"
FASTQ_R1="<sample_R1.fastq.gz>"
FASTQ_R2="<sample_R2.fastq.gz>"
TRIM_DIR="<trimmed_reads_directory>"
QC_DIR="<fastqc_output_directory>"
ALIGN_DIR="<alignment_output_directory>"
STAR_INDEX="<GRCh38.p14_Ensembl_v110_STAR_index>"
RSEQC_BED="<matching_transcriptome_annotation.bed12>"
ANNOTATION_GTF="<GRCh38.p14_Ensembl_v110_annotation.gtf>"

TRIMMED_R1="${TRIM_DIR}/${SAMPLE_ID}_R1.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_ID}_R2.fastq.gz"
BAM="${ALIGN_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam"

mkdir -p "$TRIM_DIR" "$QC_DIR" "$ALIGN_DIR"

# Trim adapters and filter low-quality reads with fastp v0.23.4 (default settings).
fastp \
  --in1 "$FASTQ_R1" --in2 "$FASTQ_R2" \
  --out1 "$TRIMMED_R1" --out2 "$TRIMMED_R2"

# Inspect trimmed-read quality with FastQC v0.11.9 (default settings).
fastqc --outdir "$QC_DIR" "$TRIMMED_R1" "$TRIMMED_R2"

# Align to the pre-built GRCh38.p14 / Ensembl v110 STAR index with STAR v2.7.10b.
# Compressed FASTQ input and coordinate-sorted BAM output are specified for downstream tools.
STAR \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$TRIMMED_R1" "$TRIMMED_R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE_ID}_" \
  --outSAMtype BAM SortedByCoordinate

# Assess alignment quality with RSeQC v5.0.1 using a matching BED12 annotation.
infer_experiment.py -r "$RSEQC_BED" -i "$BAM"
read_distribution.py -r "$RSEQC_BED" -i "$BAM"
geneBody_coverage.py -r "$RSEQC_BED" -i "$BAM" \
  -o "${ALIGN_DIR}/${SAMPLE_ID}_gene_body_coverage"

# Count reads with Rsubread v2.20.0; paired-end mode is enabled for this example.
# Other featureCounts options retain their defaults (including unstranded counting).
Rscript - "$BAM" "$ANNOTATION_GTF" "${ALIGN_DIR}/${SAMPLE_ID}_counts.txt" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
counts <- Rsubread::featureCounts(
  files = args[1],
  annot.ext = args[2],
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE
)
write.table(
  counts$counts,
  file = args[3],
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
RSCRIPT
