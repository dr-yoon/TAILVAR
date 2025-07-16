# Load necessary libraries for data manipulation and modeling
library(tidyverse)
library(caret)
library(randomForest)
library(corrplot)
library(reshape2)
library(ggbreak)

# Load pre-processed datasets for model development and validation
HGMD_data   <- read_tsv("stoploss_HGMD_preprocessed.txt") %>% mutate(Class = "P/LP")
HGMD_tango   <- read_tsv("stoploss_HGMD_TANGO_score.tsv")
HGMD_canya   <- read_tsv("stoploss_HGMD_CANYA_score.tsv")
HGMD_data <- HGMD_data %>% left_join(HGMD_tango, by = c("Var_ID" = "Variant")) %>% left_join(HGMD_canya, by = c("Var_ID" = "seqid"))

gnomad_data   <- read_tsv("stoploss_gnomad_preprocessed.txt") %>% mutate(Class = "B/LB")
gnomad_tango   <- read_tsv("stoploss_gnomad_TANGO_score.tsv")
gnomad_canya   <- read_tsv("stoploss_gnomad_CANYA_score.tsv")
gnomad_data <- gnomad_data %>% left_join(gnomad_tango, by = c("Var_ID" = "Variant")) %>% left_join(gnomad_canya, by = c("Var_ID" = "seqid"))
gnomad_data <- gnomad_data %>% filter(!Clinvar %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Conflicting_classifications_of_pathogenicity"))

train_data <- bind_rows(HGMD_data, gnomad_data) %>% distinct() # Model training set

Clinvar_data <- read_tsv("stoploss_Clinvar_preprocessed.txt")
Clinvar_data <- Clinvar_data %>% filter(Clinvar %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
                                                       "Uncertain_significance", "Benign", "Likely_benign", "Benign/Likely_benign")) %>%
  mutate(Class = case_when(str_detect(Clinvar, regex("athogenic", ignore_case = TRUE)) ~ "P/LP",
                           str_detect(Clinvar, regex("ncertain",  ignore_case = TRUE)) ~ "VUS", str_detect(Clinvar, regex("enign", ignore_case = TRUE)) ~ "B/LB", TRUE ~ Clinvar))
Clinvar_tango   <- read_tsv("stoploss_Clinvar_TANGO_score.tsv")
Clinvar_canya   <- read_tsv("stoploss_Clinvar_CANYA_score.tsv")
Clinvar_data <- Clinvar_data %>% left_join(Clinvar_tango, by = c("Var_ID" = "Variant")) %>% left_join(Clinvar_canya, by = c("Var_ID" = "seqid"))
Clinvar_data <- Clinvar_data %>% filter(Class %in% c("B/LB", "P/LP"))

RME_data   <- read_tsv("stoploss_RME_preprocessed.txt") %>% mutate(Class = "B/LB")
RME_tango   <- read_tsv("stoploss_RME_TANGO_score.tsv")
RME_canya   <- read_tsv("stoploss_RME_CANYA_score.tsv")
RME_data <- RME_data %>% left_join(RME_tango, by = c("Var_ID" = "Variant")) %>% left_join(RME_canya, by = c("Var_ID" = "seqid"))
RME_data <- RME_data %>% filter(!Clinvar %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Conflicting_classifications_of_pathogenicity"))

test_data  <- bind_rows(Clinvar_data, RME_data) %>% distinct()  # Model testing set

write.table(train_data, "Train_data_input.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(test_data, "Test_data_input.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Dataset for predicting TAILVAR score on all stoploss SNVs
SNV_data    <- read_tsv("stoploss_SNV_preprocessed.txt")
SNV_tango   <- read_tsv("stoploss_SNV_TANGO_score.tsv")
SNV_canya   <- read_tsv("stoploss_SNV_CANYA_score.tsv")
SNV_data <- SNV_data %>% left_join(SNV_tango, by = c("Var_ID" = "Variant")) %>% left_join(SNV_canya, by = c("Var_ID" = "seqid"))

# Define feature groups for model input
Numeric_variables <- c("CADD", "GERP", "phyloP100", "phastCons100", "pLI", "LOEUF", "UTR3_length", "UTR3_GC", "protein_length", "extended_length", "hydro_KD", "hydro_MJ","mRNA_stability", "TANGO", "CANYA")
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Impute missing values with the median (for numeric data)
impute_median <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df) && is.numeric(df[[v]])) {
      df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE)
    }
  }
  return(df)
}

train_data <- impute_median(train_data, Numeric_variables)
test_data <- impute_median(test_data, Numeric_variables)
SNV_data <- impute_median(SNV_data, Numeric_variables)

# Prepare the dataset
train_data <- train_data %>% filter(Class %in% c("B/LB", "P/LP")) %>% mutate(Class = recode(Class, "B/LB" = "B_LB", "P/LP" = "P_LP"), Class = factor(Class, levels = c("B_LB", "P_LP")))
test_data <- test_data %>% filter(Class %in% c("B/LB", "P/LP")) %>% mutate(Class = recode(Class, "B/LB" = "B_LB", "P/LP" = "P_LP")) %>% dplyr::select(-Upload_variation)

Model_variables  <- c("Upload_variation", "Class", all_of(Numeric_variables), all_of(amino_acid_columns), "IDP")
comp_scores <- c("CADD", "DANN", "Eigen_PC", "BayesDel_noAF", "BayesDel_addAF", "FATHMM_MKL", "fitCons_int", "MutationTaster", "VEST4", "GPN_MSA", "GERP", "phyloP100", "phastCons100")

train_data <- impute_median(train_data, comp_scores)
test_data <- impute_median(test_data, comp_scores)
SNV_data <- impute_median(SNV_data, comp_scores)

# Random Forest model development
set.seed(123)
model_train <- train_data %>% dplyr::select(all_of(Model_variables), -Upload_variation)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define grid of hyper-parameters to tune (mtry, ntree)
tune_grid <- expand.grid(mtry = c(2, 5, 10, 15, 20))
ntree_values <- c(20, 50, 100, 200, 300, 400)
maxnodes_values <- c(5, 10, 20, 30)

results <- data.frame()

# Perform grid search to optimize mtry, ntree, and maxnodes
for (maxnodes in maxnodes_values) {
  for (ntree in ntree_values) {
    # Train the model with current ntree value and the entire tune grid for mtry
    rf_tune <- train(Class ~ ., data = model_train, method = "rf", 
                     trControl = train_control, tuneGrid = tune_grid, metric = "ROC", 
                     ntree = ntree, maxnodes = maxnodes)
    
    # Loop through all combinations of mtry to extract results
    for (i in 1:nrow(tune_grid)) {
      mtry <- tune_grid$mtry[i]
      
      # Find the ROC for the current combination of mtry
      roc_value <- rf_tune$results$ROC[rf_tune$results$mtry == mtry]
      
      # Store the result
      results <- rbind(results, data.frame(ntree = ntree, mtry = mtry, maxnodes = maxnodes, ROC = roc_value))
    }
  }
}

# Select the best AUROC combination of mtry, ntree, maxnodes
selected_params <- c(20, 100, 30)
cat("Selected mtry:", selected_params[1], "ntree:", selected_params[2], "maxnode:", selected_params[3],"\n")

# Train the final Random Forest model using the selected hyper-parameters
final_rf_model <- randomForest(Class ~ ., data = model_train, mtry = selected_params[1], ntree = selected_params[2], maxnodes = selected_params[3])

# Predict TAILVAR scores for all datasets using the final model
train_data$TAILVAR <- predict(final_rf_model, train_data, type = "prob")[, "P_LP"]
test_data$TAILVAR <- predict(final_rf_model, test_data, type = "prob")[, "P_LP"]
SNV_data$TAILVAR <- predict(final_rf_model, SNV_data, type = "prob")[, "P_LP"]

write.table(train_data, "Train_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(test_data, "Test_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(SNV_data, "stoploss_SNV_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
