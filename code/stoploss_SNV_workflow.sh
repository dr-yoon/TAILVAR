#!/bin/bash
set -euo pipefail  # Exit on error, treat unset variables as an error, and pipe failures propagate

# Define working directory and script paths
WORK_DIR="/path_to/TAILVAR"
cd "${WORK_DIR}"

# Step 1: Extract stop codon positions from GENCODE
GENCODE="gencode.v47.annotation.gtf.gz" # Path to GENCODE annotation file
FILE_NAME="stoploss_SNV"
VCF_FILE="${FILE_NAME}.vcf"
HEADER="vcf_header.txt"

awk -F '\t' '
BEGIN {
  OFS = "\t"
  print "chromosome", "start", "end", "strand", "gene_id", "transcript_id", "gene_name"
}
$3 == "stop_codon" && $9 ~ /tag "MANE_Select"/ {
  match($9, /gene_id "([^"]+)"/, a)
  match($9, /transcript_id "([^"]+)"/, b)
  match($9, /gene_name "([^"]+)"/, c)
  print $1, $4, $5, $7, a[1], b[1], c[1]
}
' <(zcat "$GENCODE") > MANE_stop_codon_info.tsv

# Step 2: Generate a VCF file with all-possible SNVs at stop codons
VCF_BODY="${FILE_NAME}.tsv"
Rscript extract_transcript_info.R ${VCF_BODY}

tr -d '\r' < "${VCF_BODY}" | \
awk '{print $1, $2, ".", $3, $4, ".", ".", "."}' OFS="\t" >> ${FILE_NAME}_tmp.txt
cat "${HEADER}" "${FILE_NAME}_tmp.txt" > "${VCF_FILE}"
rm -f "${FILE_NAME}_tmp.txt"

# Step 3: Annotate the VCF file using VEP (Variant Effect Predictor)
SING_IMAGE="/path_to/ensembl-vep_latest.sif"  # Path to VEP Singularity image
VEP_DATA="/path_to/vep_data"                  # Path to VEP data directory
REF_GENOME="/path_to/Homo_sapiens_assembly38.fasta"  # Path to reference genome (hg38)
DBNSFP4_DB="${VEP_DATA}/plugins/dbNSFP4.9a_grch38.gz"     # Path to dbNSFP 4.9a database
DBNSFP5_DB="${VEP_DATA}/plugins/dbNSFP5.1a_grch38.gz"     # Path to dbNSFP 5.1a database
GPN_MSA_DB="/path_to/GPN_MSA_scores.vcf.gz" # Path to GPN_MSA score database
CLINVAR_DB="/path_to/clinvar_20250623.vcf.gz" # Path to ClinVar database
gnomAD_DB="/path_to/gnomad.joint.v4.1.sites.vcf.gz" # Path to gnomAD v4.1 database

# Annotation command using VEP
singularity run "${SING_IMAGE}" vep \
  --dir "${VEP_DATA}" \
  --dir_plugins "${VEP_DATA}/plugins" \
  --cache --offline \
  --fasta "${REF_GENOME}" \
  --format vcf --force_overwrite --no_stats --tab --show_ref_allele --hgvs --symbol \
  --biotype --canonical --numbers --mane --variant_class --protein \
  --input_file "${WORK_DIR}/${VCF_FILE}" \
  --output_file "${WORK_DIR}/${FILE_NAME}_anno.txt" \
  --plugin dbNSFP,${DBNSFP5_DB},MutationTaster_score,VEST4_score,BayesDel_noAF_score,BayesDel_addAF_score,CADD_phred,DANN_score,Eigen-phred_coding,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,RegeneronME_ALL_AF,AllofUs_POPMAX_AF,ALFA_Total_AF \
  --plugin dbNSFP,${DBNSFP4_DB},fathmm-MKL_coding_score,integrated_fitCons_score \
  --custom file=${vep_data}/plugins/GPN_MSA_scores.vcf.gz,short_name=GPN_MSA,format=vcf,type=exact,fields=score \
  --custom file=${clinvar_db},short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG \
  --custom file=${gnomad_db},short_name=gnomAD4.1,format=vcf,type=exact,fields=nhomalt_joint%AF_grpmax_joint%AF_joint

# Step 4: Filter and preprocess annotated variants
Rscript filter_variants.R "${FILE_NAME}"
Rscript TAILVAR_preprocess.R "${FILE_NAME}"

# Step 5: Run aggregation prediction tools, TANGO and CANYA
bash run_aggregation.sh "${WORK_DIR}" "${FILE_NAME}"
