# TAILVAR: Development and Assessment

## Overview

### Training set
- **`P/LP'**: HGMD 'DM' category
- **`B/LB'**: [gnomAD](https://gnomad.broadinstitute.org/) AF > 0.1% or homozygous allele

### Test set
- **`P/LP'**: ClinVar 'P/LP' category
- **`B/LB'**: AF > 0.1% in public databases ([ALFA](https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/), [Regeneron Million Exome](https://rgc-research.regeneron.com/), [AllofUs Genome](https://support.researchallofus.org/))

#### Input Files
- **`Preprocessed data`**: Preprocessed annotation files.
- **`TANGO and CANYA predictions`**: Aggregation predictions for C-terminal peptides.

### Random Forest Model
#### Script: `TAILVAR_score.R`
This R script builds a Random Forest model to predict the pathogenicity of stoploss variants and calculate the TAILVAR scores.

#### Required Packages
- `tidyverse`, `caret`, `randomForest`, `readr`

#### Files
- **`Train_data_input.txt`**: Training dataset for model development.
- **`Test_data_input.txt`**: Testing dataset for model validation.

#### Output Files
- **`Train_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model development dataset.
- **`Test_TAILVAR_score.txt`**: Calculated TAILVAR scores for the model validation dataset.
- **`stoploss_SNV_TAILVAR_score.txt`**: Calculated TAILVAR scores for all possible stoploss variants.
