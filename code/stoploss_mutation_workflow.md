# All-Possible Single-Nucleotide Mutations at Stop Codons

This repository provides a pipeline for extracting, annotating, and filtering variants at stop codons, focusing on all possible single-nucleotide substitutions (SNVs) or insertion/deletion. The process involves querying the Ensembl database, generating a VCF file, annotating it using the Variant Effect Predictor (VEP), and filtering the results.

## Pipeline Overview

### Step 1: Extract Stop Codon Information from the **[GENCODE] https://www.gencodegenes.org/)** 
#### Output Files:
- **`MANE_stop_codon_info.tsv`**: Information on stop codons for all MANE select transcripts retrieved from the GENCODE database.

### Step 2: Generate a VCF file with all-possible mutations at stop codons
#### Script: `stop_codon_positions.R`
This R script extracts information on stop codon positions from MANE transcripts.

#### Required packages:
To run this script, you need to have the following R packages installed:
- `tidyverse`

This step will generate a vcf file for all possible single-nucleotide substitutions (SNV) or insertions/deletions (INS/DEL).
#### Output Files:
- **`stoploss_SNV.vcf`**: a vcf file with all possible single-nucleotide substitutions (SNV).
- **`stoploss_DEL.vcf`**: a vcf file with all possible single-nucleotide deletions (DEL).
- **`stoploss_INS.vcf`**: a vcf file with all possible single-nucleotide insertions (INS).

### Step 3: Annotate the VCF File Using VEP (Variant Effect Predictor)
This step will annotate the vcf files above to add additional information about the variants.

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

#### Output Files: Annotated VCF file with variant information.
- **`stoploss_SNV_anno.txt`**
- **`stoploss_DEL_anno.txt`**
- **`stoploss_INS_anno.txt`**

### Step 4: Filter and preprocess annotated variants
#### Script: `filter_variants.R`
This R script filters out unnecessary columns and formats the annotated file to retain only relevant information.

#### Output Files: The final filtered and formatted variant list.
- **`stoploss_SNV_filtered.txt`**
- **`stoploss_DEL_filtered.txt`**
- **`stoploss_INS_filtered.txt`**

#### Script: `TAILVAR_preprocess.R`
This R script will add required information (protein information on translated C-terminal extensions) and formats the file to retain only required field for machine learning.

#### Required packages:
To run this script, you need to have the following R packages installed:
- `tidyverse`,`biomaRt`,`stringr`,`Biostrings`,`Peptides`

#### Output Files: The preprocessed input file for machine learning model.
- **`stoploss_SNV_preprocessed.txt`**
- **`stoploss_DEL_preprocessed.txt`**
- **`stoploss_INS_preprocessed.txt`**

### Step 5: Predict aggregation properties using TANGO and CANYA
#### Script: `run_aggregation.sh`
This script execute to run TANGO and CANYA for C-terminal peptides to predict the aggregation properties.

Academic License should be obtained for [TANGO](https://tango.crg.es/)
Conda environment should be install for [CANYA](https://github.com/lehner-lab/canya)


## How to Run the Pipeline

1. Clone this repository or download the scripts.
2. Ensure all required R packages, tools, and data files are installed on your system.
3. Update the paths in the Bash script (`stoploss_SNV_workflow.sh`) to match your directory structure.
4. Run the pipeline by executing the Bash script:

   ```bash
   ./stoploss_mutation_workflow.sh

