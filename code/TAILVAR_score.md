# TAILVAR Score: Model Development and Validation

This repository provides an overview of the scripts used for the development and validation of the TAILVAR score, a tool for predicting the pathogenicity of stoplost variants.

## Overview

### Step 1: Prepare Input Files
#### Script: `prepare_TAILVAR_input.R`
This R script is responsible for combining and formatting data necessary for downstream analysis in the TAILVAR score pipeline.

#### Required Packages
- `tidyverse`

#### Input Files
- **`mane_transcripts_info.txt`**: Contains information on stop codons for all MANE select transcripts retrieved from the Ensembl database.
- **`stoplost_SNV_filtered.txt`**: An annotated and filtered file containing all possible single-nucleotide substitutions at stop codons.
- **`stoplost_HGMD_filtered.txt`**: An annotated and filtered file of HGMD stoplost variants.
- **`stoplost_ClinVar_filtered.txt`**: An annotated and filtered file of ClinVar stoplost variants.

#### Output Files
- **`HGMD_gnomAD_train_data.txt`**: Training dataset for model development.
- **`ClinVar_ALFA_test_data.txt`**: Testing dataset for model validation.
- **`stoplost_SNV_prediction_data.txt`**: Input data for TAILVAR score prediction of all possible stoplost variants.

### Step 2: Build Ensemble Model
#### Script: `TAILVAR_score.R`
This R script builds a Random Forest model to predict the pathogenicity of stoplost variants and calculate the TAILVAR scores.

#### Required Packages
- `tidyverse`, `caret`, `randomForest`, `readr`, `corrplot`, `reshape2`, `ggplot2`

#### Input Files
- **`HGMD_gnomAD_train_data.txt`**: Training dataset for model development.
- **`ClinVar_ALFA_test_data.txt`**: Testing dataset for model validation.
- **`stoplost_SNV_prediction_data.txt`**: Input data for TAILVAR score prediction of all possible stoplost variants.

#### Output Files
- **`Training_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model development dataset.
- **`Testing_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model validation dataset.
- **`stoplost_SNV_TAILVAR_score.txt`**: Calculated TAILVAR scores for all possible stoplost variants.
- **`Correlation_plot.svg`**
- **`Parameter_Importance.svg`**

### Step 3: Generate Performance Plots
#### Script: `TAILVAR_performance.R`
This R script generates plots for performance comparisons with TAILVAR score.

#### Required Packages
- `tidyverse`, `pROC`, `ggplot2`, `svglite`

#### Input Files
- **`Training_TAILVAR_score.txt`**: Calculated TAILVAR scores for training dataset.
- **`Testing_TAILVAR_score.txt`**: Calculated TAILVAR scores for testing dataset.
- **`stoplost_SNV_TAILVAR_score.txt`**: Calculated TAILVAR scores for all possible stoplost variants.

#### Output Files
- **`Development_dataset_TAILVAR_distribtion.svg`**
- **`Validation_dataset_TAILVAR_distribtion.svg`**
- **`gnomAD_AF_TAILVAR_relationship.svg`**
- **`Development_dataset_TAILVAR_AUROC.svg`**
- **`Validation_dataset_TAILVAR_AUROC.svg`**
- **`TAILVAR_AUROC_comparisons.svg`**

