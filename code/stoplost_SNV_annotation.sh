#!/bin/bash
set -euo pipefail  # Exit on error, treat unset variables as an error, and pipe failures propagate

# Define working directory and script paths
WORK_DIR="/path_to/TAILVAR/code"
cd "${WORK_DIR}"

# Step 1: Extract variants and information from Ensembl database
Rscript extract_ensembl_data.R

# Step 2: Generate a VCF file with saturated mutations at stop codons
FILE_NAME="stoplost_SNV"
MUTA_SITES="${FILE_NAME}.txt"
VCF_FILE="${FILE_NAME}.vcf"
HEADER="vcf_header.txt"

# Convert to Unix line endings, create VCF body, and concatenate with header
tr -d '\r' < "${MUTA_SITES}" | \
awk '{print $1, $2, ".", $3, $4, ".", ".", "."}' OFS="\t" >> "${FILE_NAME}_body.txt"

cat "${HEADER}" "${FILE_NAME}_body.txt" > "${VCF_FILE}"
rm -f "${FILE_NAME}_body.txt"  # Clean up the temporary file

# Step 3: Annotate the VCF file using VEP (Variant Effect Predictor)
SING_IMAGE="/path_to/ensembl-vep_latest.sif"  # Path to VEP Singularity image
VEP_DATA="/path_to/vep_data"                  # Path to VEP data directory
REF_GENOME="/path_to/Homo_sapiens_assembly38.fasta"  # Path to reference genome (hg38)
DBNSFP_DB="${VEP_DATA}/plugins/dbNSFP4.9a_grch38.gz"     # Path to dbNSFP database
CLINVAR_DB="/path_to/clinvar_20241021.vcf.gz" # Path to ClinVar database

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
  --plugin dbNSFP,"${DBNSFP_DB}",CADD_phred,DANN_score,fathmm-MKL_coding_score,Eigen-phred_coding,BayesDel_addAF_score,BayesDel_noAF_score,integrated_fitCons_score,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,gnomAD_exomes_flag,gnomAD_exomes_AF,gnomAD_exomes_POPMAX_AF,gnomAD_genomes_flag,gnomAD_genomes_AF,gnomAD_genomes_POPMAX_AF, ALFA_European_AF,ALFA_African_AF,ALFA_Asian_AF,ALFA_Other_AF,ALFA_Total_AF \
  --custom file="${CLINVAR_DB}",short_name=ClinVar,format=vcf,type=exact,fields=CLNDN%CLNSIG

# Step 4: Filter annotated variants
Rscript filter_variants.R "${FILE_NAME}"
