# All-Possible Single-Nucleotide Substitutions (SNVs) at Stop Codons

This repository provides a pipeline for extracting, annotating, and filtering variants at stop codons, focusing on all possible single-nucleotide substitutions (SNVs). The process involves querying the Ensembl database, generating a VCF file, annotating it using the Variant Effect Predictor (VEP), and filtering the results.

## Pipeline Overview

### Step 1: Extract Variants and Information from the Ensembl Database
#### Script: `extract_ensembl_data.R`
This R script extracts data from the Ensembl database, focusing on stop-lost variants and associated transcript information.

#### Required packages:
To run this script, you need to have the following R packages installed:
- `tidyverse`, `biomaRt`, `stringr`, `Biostrings`

#### Output Files:
- **`mane_transcripts_info.txt`**: Information on stop codons for all MANE select transcripts retrieved from the Ensembl database.
- **`stoplost_SNV.txt`**: A list of all possible single-nucleotide substitutions at stop codons, including chromosome, position, reference, and alternative sequences.

### Step 2: Generate a VCF File
The **`stoplost_SNV.txt`** file generated in Step 1 is converted to VCF format, resulting in the **`stoplost_SNV.vcf`** file.

### Step 3: Annotate the VCF File Using VEP (Variant Effect Predictor)
The **`stoplost_SNV.vcf`** file is annotated using VEP to add additional information about the variants.

#### Requirements:
To run the VEP annotation, you need the following tools and data:
- **[Singularity](https://github.com/sylabs/singularity)**: Used to run VEP in a containerized environment. Alternatively, VEP can be run directly without Singularity.
- **[VEP (Variant Effect Predictor)](https://github.com/Ensembl/ensembl-vep)**
- **VEP Plugins and Database Files**:
  - [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)
  - [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
  - [NARD2](https://nard.macrogen.com/)
  - [ToMMo 60KJPN](https://jmorp.megabank.tohoku.ac.jp/downloads/tommo-60kjpn-20240904-af_snvindelall)
  - Reference genome file (hg38)

#### Output Files:
- **`stoplost_SNV_anno.txt`**: Annotated VCF file with variant information.

### Step 4: Filter Annotated Variants
#### Script: `filter_variants.R`
This R script filters out unnecessary columns and formats the annotated file to retain only relevant information.

#### Output Files:
- **`stoplost_SNV_filtered.txt`**: The final filtered and formatted variant list.

## How to Run the Pipeline

1. Clone this repository or download the scripts.
2. Ensure all required R packages, tools, and data files are installed on your system.
3. Update the paths in the Bash script (`stoplost_SNV_annotation.sh`) to match your directory structure.
4. Run the pipeline by executing the Bash script:

   ```bash
   ./stoplost_SNV_annotation.sh

#### Output Summary
- **`mane_transcripts_info.txt`**: Stop codon information for MANE select transcripts.
- **`stoplost_SNV.txt`**: All possible SNVs at stop codons.
- **`stoplost_SNV.vcf`**: VCF file containing SNVs.
- **`stoplost_SNV_anno.txt`**: Annotated VCF file.
- **`stoplost_SNV_filtered.txt`**: Final filtered variant list.
