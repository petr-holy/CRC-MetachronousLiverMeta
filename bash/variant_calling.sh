# Reporting summary only; placeholders describe input files and are not study-specific.

REFERENCE=reference_GRCh38.fa
TARGETS=target_regions.interval_list
TUMOR_BAM=tumor_sample.bam
NORMAL_BAM=normal_sample.bam
TUMOR_ID=tumor_sample
NORMAL_ID=normal_sample
GNOMAD_AF=gnomad_population_AF.vcf.gz
PON=somatic_panel_of_normals.vcf.gz
DBSNP=dbSNP151.vcf.gz

# Somatic contamination and artifact estimation (GATK).
gatk GetPileupSummaries \
  --reference "$REFERENCE" --intervals "$TARGETS" --input "$TUMOR_BAM" \
  --variant "$GNOMAD_AF" --interval-set-rule INTERSECTION \
  --output tumor_pileups.table

gatk CalculateContamination \
  --input tumor_pileups.table --tumor-segmentation tumor_segments.txt \
  --output tumor_contamination.table

gatk CollectSequencingArtifactMetrics \
  --REFERENCE_SEQUENCE "$REFERENCE" --INPUT "$TUMOR_BAM" \
  --OUTPUT tumor_artifact_metrics --FILE_EXTENSION .txt

# Somatic SNVs and indels: Mutect2, VarDict, VarScan, and Strelka.
gatk Mutect2 \
  --reference "$REFERENCE" --intervals "$TARGETS" \
  --input "$TUMOR_BAM" --input "$NORMAL_BAM" -normal "$NORMAL_ID" \
  --germline-resource "$GNOMAD_AF" --panel-of-normals "$PON" \
  --f1r2-tar-gz tumor_f1r2.tar.gz --native-pair-hmm-threads 60 \
  --output mutect2_raw.vcf.gz

gatk LearnReadOrientationModel \
  --input tumor_f1r2.tar.gz --output orientation_model.tar.gz

gatk FilterMutectCalls \
  --reference "$REFERENCE" --variant mutect2_raw.vcf.gz \
  --tumor-segmentation tumor_segments.txt \
  --contamination-table tumor_contamination.table \
  --orientation-bias-artifact-priors orientation_model.tar.gz \
  --output mutect2_filtered.vcf.gz


bcftools filter -i 'STATUS ~"StrongSomatic"' vardict_paired.vcf \
  -Oz -o vardict_strong_somatic.vcf.gz
bcftools filter -i 'MQ>=30 & FORMAT/DP>=30 & FORMAT/VD>=5 & QUAL>=100' \
  vardict_strong_somatic.vcf.gz -Oz -o vardict_filtered.vcf.gz

samtools mpileup --verbosity 5 -q 1 -B -f "$REFERENCE" \
  -o tumor_normal.pileup "$NORMAL_BAM" "$TUMOR_BAM"
varscan somatic tumor_normal.pileup varscan_output \
  --mpileup 1 --min-coverage 8 --min-var-freq 0.05 \
  --min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.0 \
  --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf
varscan processSomatic varscan_output.snp.vcf \
  --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07
varscan processSomatic varscan_output.indel.vcf \
  --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07

configManta.py --normalBam "$NORMAL_BAM" --tumorBam "$TUMOR_BAM" \
  --referenceFasta "$REFERENCE" --exome --runDir manta_run
manta_run/runWorkflow.py --jobs 60 --memGb 100
configureStrelkaSomaticWorkflow.py \
  --normalBam "$NORMAL_BAM" --tumorBam "$TUMOR_BAM" \
  --referenceFasta "$REFERENCE" --indelCandidates manta_candidate_indels.vcf.gz \
  --runDir strelka_run
strelka_run/runWorkflow.py --mode local --jobs 60 --memGb 100

# Normalize called VCFs and combine caller evidence.
vt normalize caller_variants.vcf.gz -r "$REFERENCE" -o normalized_variants.vcf
java -jar somaticCombiner.jar --threshold 0.05 \
  --mutect2 mutect2_filtered_normalized.vcf \
  --varscan-snv varscan_snv_normalized.vcf \
  --varscan-indel varscan_indel_normalized.vcf \
  --strelka-snv strelka_snvs_normalized.vcf \
  --strelka-indel strelka_indels_normalized.vcf \
  -o combined_somatic.vcf
bcftools filter -i 'FILTER ~"ADJ_PASS"' combined_somatic.vcf \
  -o combined_somatic_ADJ_PASS.vcf

# Somatic annotation and VCF-to-MAF conversion (VEP cache 111; GRCh38).
vcf2maf.pl --input-vcf combined_somatic_ADJ_PASS.vcf \
  --output-maf somatic_variants.maf --tumor-id "$TUMOR_ID" \
  --normal-id "$NORMAL_ID" --vep-path vep_executable --vep-data vep_cache \
  --max-subpop-af 0.001 --ncbi-build GRCh38 --ref-fasta "$REFERENCE" \
  --vep-forks 4 --cache-version 111 --vep-overwrite --verbose \
  --vcf-tumor-id "$TUMOR_ID" --vcf-normal-id "$NORMAL_ID" \
  --retain-info NumCallers,Tumor_AF,Tumor_DP

