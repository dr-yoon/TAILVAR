# Load necessary libraries for data manipulation and modeling
library(tidyverse)
library(caret)
library(randomForest)
library(corrplot)
library(reshape2)
library(ggbreak)

# Load pre-processed datasets for model development and validation
train_data <- read_tsv("Train_data_input.txt")        # Model training set
test_data <- read_tsv("Test_data_input.txt")  # Model testing set
stoploss_all <- read_tsv("stoploss_SNV_input.txt") # Dataset for predicting TAILVAR score on all stoploss variants

# Define feature groups for model input
comp_scores <- c("CADD", "DANN", "FATHMM", "Eigen", "BayesDel_addAF", "BayesDel_noAF", "integrated_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_GC", "UTR3_length", "Extension_lengths")
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Hydrophobicity")
amino_acid_columns_rename <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "H_index")

train_data <- train_data %>% rename_at(vars(amino_acid_columns), ~ amino_acid_columns_rename)
test_data <- test_data %>% rename_at(vars(amino_acid_columns), ~ amino_acid_columns_rename)
stoploss_all <- stoploss_all %>% rename_at(vars(amino_acid_columns), ~ amino_acid_columns_rename)

# Prepare the input data for Random Forest model development
train_data <- train_data %>% 
  dplyr::select("Class", all_of(comp_scores), all_of(gene_feature), all_of(amino_acid_columns_rename)) %>% 
  mutate(Class = factor(Class, levels = c("B", "P")))  # Ensure the target variable is a factor with correct levels

# Set up cross-validation for model training
set.seed(123)  # Set seed for reproducibility
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
options(digits = 10)

# Define grid of hyper-parameters to tune (mtry, ntree)
tune_grid <- expand.grid(mtry = c(2, 5, 10, 15, 20))  # Different values of mtry (number of variables considered at each split)
ntree_values <- c(20, 40, 60, 80, 100)                # Different values of ntree (number of trees in the forest)
maxnodes_values <- c(5, 10, 20, 30)                      # Values for maxnodes

# Initialize a dataframe to store tuning results
results <- data.frame()

# Perform grid search to optimize mtry, ntree, and maxnodes
for (maxnodes in maxnodes_values) {
  for (ntree in ntree_values) {
    # Train the model with current ntree value and the entire tune grid for mtry
    rf_tune <- train(Class ~ ., data = train_data, method = "rf", 
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

# Plot AUROC vs. hyper-parameters (mtry, ntree, maxnodes) to visualize model performance
svg("Hyperparameters.svg", width = 8, height = 8)
ggplot(results, aes(x = ntree, y = ROC, color = factor(mtry))) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "AUROC vs. Hyper-parameters (ntree, mtry) Faceted by maxnodes",
       x = "Number of Trees (ntree)", y = "AUROC", color = "mtry") +
  facet_wrap(~ maxnodes) +  # Add faceting by maxnodes
  theme_minimal() + theme(legend.position = "right")
dev.off()

# Select the best AUROC combination of mtry, ntree, maxnodes
selected_params <- c(10, 100, 30)
cat("Selected mtry:", selected_params[1], "ntree:", selected_params[2], "maxnode:", selected_params[3],"\n")

# Train the final Random Forest model using the selected hyper-parameters
final_rf_model <- randomForest(Class ~ ., data = train_data, 
                               mtry = selected_params[1], ntree = selected_params[2], maxnodes = selected_params[3])

# Predict TAILVAR scores for all datasets using the final model
train_data$TAILVAR <- predict(final_rf_model, train_data, type = "prob")[, "P"]
test_data$TAILVAR <- predict(final_rf_model, test_data, type = "prob")[, "P"]
stoploss_all$TAILVAR <- predict(final_rf_model, stoploss_all, type = "prob")[, "P"]
stoploss_all_vcf <- stoploss_all %>% select(Chromosome, Position, REF_ALLELE, Allele, TAILVAR)

# Save the TAILVAR scores to output files
write.table(train_data, "Train_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(test_data, "Test_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoploss_all, "stoploss_SNV_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoploss_all_vcf, "TAILVAR_score_vcf_input.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# Calculate the correlation matrix for the selected parameters
correlation_matrix <- cor(train_data %>% dplyr::select(all_of(comp_scores),all_of(gene_feature),all_of(gene_feature),"H_index"), method = "spearman", use = "complete.obs")
print(correlation_matrix)
correlation_melt <- melt(correlation_matrix)

# Plot the Spearman rank correlation matrix using corrplot
col <- colorRampPalette(c("#4477AA","#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
svg("Correlation_plot.svg", width = 8, height = 8)
corrplot(correlation_matrix, method = "circle", type = "upper", col = col(100), 
         tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, mar = c(0, 0, 1, 0))
dev.off()

# Extract variable importance using Gini index
variable_importance <- importance(final_rf_model, type = 2)
importance_df <- data.frame(Variable = rownames(variable_importance), Importance = variable_importance[, "MeanDecreaseGini"])
importance_df <- importance_df %>% mutate(RelativeImportance = Importance / sum(Importance) * 100) %>%
  arrange(desc(RelativeImportance)) %>%  slice_head(n = 20)

# Plot the normalized importance scores
svg("Feature_Importance.svg", width = 9, height = 6)
importance_df %>% arrange(desc(RelativeImportance)) %>%
  ggplot(aes(x = reorder(Variable, RelativeImportance, decreasing = TRUE), y = RelativeImportance)) +
  geom_bar(stat = "identity", fill = "#228B22") + geom_text(aes(label = round(RelativeImportance, 1)), vjust = 1.5, size = 3.5, color = "black") +
  labs(x = "", y = "Relative importance (%)") +
  theme_classic() + theme(axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 14, angle = -45, hjust = 0)) + scale_y_break(c(10, 30))
dev.off()
