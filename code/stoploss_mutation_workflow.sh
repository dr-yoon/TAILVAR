#!/bin/bash
set -euo pipefail  # Exit on error, treat unset variables as an error, and pipe failures propagate

# Define working directory and script paths
WORK_DIR="/path_to/TAILVAR"
cd "${WORK_DIR}"

# Step 1: Extract stop codon positions from GENCODE
GENCODE="gencode.v47.annotation.gtf.gz" # Path to GENCODE annotation file
FILE_NAME="stoploss_SNV" # or "stoploss_DEL" / "stoploss_INS", change this by variant_type
VCF_FILE="${FILE_NAME}.vcf"
HEADER="vcf_header.txt"
var_type="SNV" # or "DEL" / "INS", change this by variant_type

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
Rscript stop_codon_positions.R ${file_name}

awk -F'\t' '{print $1":"$2"-"$2"\t"$0}' ${file_name}.tsv | \
while read region chr pos; do
    base=$(samtools faidx ${REF} $region | tail -n +2 | tr -d '\n')
    echo -e "${chr}\t${pos}\t${base}"
done > ${file_name}_ref.tsv

case "$var_type" in
    "SNV"|"DEL"|"INS")
        ;;
    *)
        echo "Warning: unknown var_type '$var_type'. Defaulting to SNV."
        var_type="SNV"
        ;;
esac

> ${file_name}_body.txt
if [[ "$var_type" == "SNV" ]]; then
    while IFS=$'\t' read -r chrom pos ref; do
        for alt in A C G T; do
            if [[ "$alt" != "$ref" ]]; then
                echo -e "${chrom}\t${pos}\t.\t${ref}\t${alt}\t.\t.\t." >> ${file_name}_body.txt
            fi
        done
    done < ${file_name}_ref.tsv

elif [[ "$var_type" == "DEL" ]]; then
    while IFS=$'\t' read -r chrom pos ref; do
        next_pos=$((pos + 1))
        region2="${chrom}:${next_pos}-${next_pos}"
        base2=$(samtools faidx ${REF} $region2 | tail -n +2 | tr -d '\n')
        if [[ "$ref" =~ ^[ACGT]$ && "$base2" =~ ^[ACGT]$ ]]; then
            echo -e "${chrom}\t${pos}\t.\t${ref}${base2}\t${base2}\t.\t.\t." >> ${file_name}_body.txt
        fi
    done < ${file_name}_ref.tsv

elif [[ "$var_type" == "INS" ]]; then
    while IFS=$'\t' read -r chrom pos ref; do
        for ins in A C G T; do
            echo -e "${chrom}\t${pos}\t.\t${ref}\t${ins}${ref}\t.\t.\t." >> ${file_name}_body.txt
        done
    done < ${file_name}_ref.tsv
fi

cat $header ${file_name}_body.txt > $vcf_file
rm -rf ${file_name}.tsv ${file_name}_ref.tsv ${file_name}_body.txt

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
