# TAILVAR: Model Development and Assessment

## Overview

### Training set
#### Script: `prepare_TAILVAR_input.R`
This R script is responsible for combining and formatting data necessary for downstream analysis in the TAILVAR score pipeline.

#### Required Packages
- `tidyverse`

#### Input Files
- **`mane_transcripts_info.txt`**: Contains information on stop codons for all MANE select transcripts retrieved from the Ensembl database.
- **`stoploss_SNV_filtered.txt`**: An annotated and filtered file containing all possible single-nucleotide substitutions at stop codons.
- **`stoploss_HGMD_filtered.txt`**: An annotated and filtered file of HGMD stoploss variants.
- **`stoploss_ClinVar_filtered.txt`**: An annotated and filtered file of ClinVar stoploss variants.

#### Output Files
- **`Train_data_input.txt`**: Training dataset for model development.
- **`Test_data_input.txt`**: Testing dataset for model validation.
- **`Stoploss_SNV_input.txt`**: Input data for TAILVAR score prediction of all possible stoploss variants.

### Step 2: Build Ensemble Model
#### Script: `TAILVAR_score.R`
This R script builds a Random Forest model to predict the pathogenicity of stoploss variants and calculate the TAILVAR scores.

#### Required Packages
- `tidyverse`, `caret`, `randomForest`, `readr`, `corrplot`, `reshape2`, `ggplot2`

#### Input Files
- **`Train_data_input.txt`**: Training dataset for model development.
- **`Test_data_input.txt`**: Testing dataset for model validation.
- **`Stoploss_SNV_input.txt`**: Input data for TAILVAR score prediction of all possible stoploss variants.

#### Output Files
- **`Train_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model development dataset.
- **`Test_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model validation dataset.
- **`stoploss_SNV_TAILVAR_score.txt`**: Calculated TAILVAR scores for all possible stoploss variants.
- **`Hyperparameters.svg`**
- **`Correlation_plot.svg`**
- **`Feature_Importance.svg`**

### Step 3: Generate Performance Plots
#### Script: `TAILVAR_performance.R`
This R script generates plots for performance comparisons with TAILVAR score.

#### Required Packages
- `tidyverse`, `pROC`, `ggplot2`, `svglite`,`mclust`

#### Input Files
- **`Train_TAILVAR_score.txt`**: Calculated TAILVAR scores for training dataset.
- **`Test_TAILVAR_score.txt`**: Calculated TAILVAR scores for testing dataset.
- **`stoploss_SNV_TAILVAR_score.txt`**: Calculated TAILVAR scores for all possible stoploss variants.

#### Output Files
- **`Train_TAILVAR_distribtion.svg`**
- **`Test_TAILVAR_distribtion.svg`**
- **`stoploss_all_GMM.svg`**
- **`Clinvar_VUS_GMM.svg`**
- **`gnomAD_AF_TAILVAR_relationship.svg`**
- **`Train_dataset_AUROC.svg`**
- **`Test_dataset_AUROC.svg`**
- **`TAILVAR_AUROC_comparisons.svg`**
