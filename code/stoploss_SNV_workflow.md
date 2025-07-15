# All-Possible Single-Nucleotide Substitutions (SNVs) at Stop Codons

This repository provides a pipeline for extracting, annotating, and filtering variants at stop codons, focusing on all possible single-nucleotide substitutions (SNVs). The process involves querying the Ensembl database, generating a VCF file, annotating it using the Variant Effect Predictor (VEP), and filtering the results.

## Pipeline Overview

### Step 1: Extract Stop Codon Information from the **[GENCODE] https://www.gencodegenes.org/)** 

### Step 2: Generate a VCF file with all-possible SNVs at stop codons
#### Script: `extract_transcript_info.R`
This R script extracts data from the Ensembl database, converting stop codon positions to a vcf file and translate 3'UTR sequences into the C-terminal extensions.

#### Required packages:
To run this script, you need to have the following R packages installed:
- `tidyverse`, `biomaRt`, `stringr`, `Biostrings`

#### Output Files:
- **`mane_transcripts_info.txt`**: Information on stop codons for all MANE select transcripts retrieved from the Ensembl database.
- **`stoploss_SNV.tsv`**: A list of all possible single-nucleotide substitutions at stop codons, including chromosome, position, reference, and alternative sequences.

### Step 3: Annotate the VCF File Using VEP (Variant Effect Predictor)
The **`stoploss_SNV.vcf`** file is annotated using VEP to add additional information about the variants.

#### Requirements:
To run the VEP annotation, you need the following tools and data files:
- **[Singularity](https://github.com/sylabs/singularity)**: Used to run VEP in a containerized environment. Alternatively, VEP can be run directly without Singularity.
- **[VEP (Variant Effect Predictor)](https://github.com/Ensembl/ensembl-vep)**
- **VEP Plugins and Database Files**: should be downloaded and placed appropriately
  - [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)
  - [GPN_MSA](https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores)
  - [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
  - [gnomAD](https://gnomad.broadinstitute.org/data#v4)
  - Reference genome file (hg38)

#### Output Files:
- **`stoploss_SNV_anno.txt`**: Annotated VCF file with variant information.

### Step 4: Filter and preprocess annotated variants
#### Script: `filter_variants.R`
This R script filters out unnecessary columns and formats the annotated file to retain only relevant information.

#### Output Files:
- **`stoploss_SNV_filtered.txt`**: The final filtered and formatted variant list.

#### Script: `TAILVAR_preprocess.R`
This R script further add required information and formats the file to retain only required field for machine learning.

#### Output Files:
- **`stoploss_SNV_preprocessed.txt`**: The preprocessed input file for machine learning model.

### Step 5: Predict aggregation properties using TANGO and CANYA
#### Script: `run_aggregation.sh`
This script execute to run TANGO and CANYA for C-terminal peptides to predict the aggregation properties.

Academic License should be obtained for [TANGO](https://tango.crg.es/)
Conda environment should be install for [CANYA](https://github.com/lehner-lab/canya)


## How to Run the Pipeline

1. Clone this repository or download the scripts.
2. Ensure all required R packages, tools, and data files are installed on your system.
3. Update the paths in the Bash script (`stoploss_SNV_annotation.sh`) to match your directory structure.
4. Run the pipeline by executing the Bash script:

   ```bash
   ./stoploss_SNV_annotation.sh

#### Output Summary
- **`mane_transcripts_info.txt`**: Stop codon information for MANE select transcripts.
- **`stoploss_SNV.txt`**: All possible SNVs at stop codons.
- **`stoploss_SNV.vcf`**: VCF file containing SNVs.
- **`stoploss_SNV_anno.txt`**: Annotated VCF file.
- **`stoploss_SNV_filtered.txt`**: Final filtered variant list.
