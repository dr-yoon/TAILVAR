# Load necessary libraries for data manipulation and modeling
library(tidyverse)
library(caret)
library(randomForest)
library(readr)
library(corrplot)
library(reshape2)

# Load pre-processed datasets for model development and validation
model_data <- read_tsv("HGMD_gnomAD_model_data.txt")        # Model development dataset
validation_data <- read_tsv("ClinVar_validation_data.txt")  # Model validation dataset
stoplost_all <- read_tsv("stoplost_SNV_prediction_data.txt")# Dataset for predicting TAILVAR score on all stoplost variants
stoplost_gnomAD <- read_tsv("stoplost_gnomAD_prediction_data.txt") # Dataset for predicting TAILVAR score on gnomAD stoplost variants

# Define feature groups for model input
comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", 
                 "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_GC", "UTR3_length", "TailAA_counts")
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Prepare the input data for Random Forest model development
model_input <- model_data %>% 
  dplyr::select("Class", all_of(comp_scores), all_of(gene_feature), all_of(amino_acid_columns)) %>% 
  mutate(Class = factor(Class, levels = c("B", "P")))  # Ensure the target variable is a factor with correct levels

# Split the data into training and testing sets (80% training, 20% testing)
set.seed(123)  # Set seed for reproducibility
trainIndex <- createDataPartition(model_input$Class, p = 0.8, list = FALSE)
train_data <- model_input[trainIndex,]
test_data <- model_input[-trainIndex,]

# Set up cross-validation for model training
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define grid of hyper-parameters to tune (mtry, ntree)
tune_grid <- expand.grid(mtry = c(2, 5, 10, 20, 30))  # Different values of mtry (number of variables considered at each split)
ntree_values <- c(10, 25, 50, 100, 150, 200)          # Different values of ntree (number of trees in the forest)

# Initialize a dataframe to store tuning results
results <- data.frame()

# Perform grid search to optimize mtry and ntree
for (ntree in ntree_values) {
  rf_tune <- train(Class ~ ., data = train_data, method = "rf", 
                   trControl = train_control, tuneGrid = tune_grid, metric = "ROC", 
                   ntree = ntree)
  
  for (mtry in tune_grid$mtry) {
    roc_value <- rf_tune$results$ROC[rf_tune$results$mtry == mtry]
    results <- rbind(results, data.frame(ntree = ntree, mtry = mtry, ROC = roc_value))
  }
}

# Plot AUROC vs. hyper-parameters (ntree, mtry) to visualize model performance
ggplot(results, aes(x = ntree, y = ROC, color = factor(mtry))) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "AUROC vs. Hyper-parameters (ntree, mtry)",
       x = "Number of Trees (ntree)", y = "AUROC", color = "mtry") +
  theme_minimal() + theme(legend.position = "right")

# Select the best combination of mtry and ntree based on the results
selected_params <- c(10, 100)  # Example selection; update based on your results
cat("Selected mtry:", selected_params[1], "with ntree:", selected_params[2], "\n")

# Train the final Random Forest model using the selected hyper-parameters
final_rf_model <- randomForest(Class ~ ., data = train_data, 
                               mtry = selected_params[1], ntree = selected_params[2])

# Predict TAILVAR scores for all datasets using the final model
model_data$TAILVAR <- predict(final_rf_model, model_data, type = "prob")[, "P"]
validation_data$TAILVAR <- predict(final_rf_model, validation_data, type = "prob")[, "P"]
stoplost_all$TAILVAR <- predict(final_rf_model, stoplost_all, type = "prob")[, "P"]
stoplost_gnomAD$TAILVAR <- predict(final_rf_model, stoplost_gnomAD, type = "prob")[, "P"]

# Save the TAILVAR scores to output files
write.table(model_data, "Development_dataset_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(validation_data, "Validation_dataset_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_gnomAD, "stoplost_gnomAD_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_all, "stoplost_SNV_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Calculate the correlation matrix for the selected parameters
correlation_matrix <- cor(model_data %>% dplyr::select(all_of(comp_scores),all_of(gene_feature)), method = "spearman", use = "complete.obs")
print(correlation_matrix)
correlation_melt <- melt(correlation_matrix)

# Plot the Spearman rank correlation matrix using corrplot
col <- colorRampPalette(c("#4477AA","#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
corrplot(correlation_matrix, method = "circle", type = "upper", col=col(100), tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, mar = c(0, 0, 1, 0))

# Extract variable importance using Gini index
variable_importance <- importance(final_rf_model, type = 2)
importance_df <- data.frame(Variable = rownames(variable_importance), Importance = variable_importance[, "MeanDecreaseGini"])
importance_df <- importance_df %>% mutate(RelativeImportance = Importance / sum(Importance) * 100) %>%
  arrange(desc(RelativeImportance)) %>%  slice_head(n = 20)

# Plot the normalized importance scores
importance_df %>% arrange(desc(RelativeImportance)) %>%
  ggplot(aes(x = reorder(Variable, RelativeImportance, decreasing = TRUE), y = RelativeImportance)) +
  geom_bar(stat = "identity", fill = "navy") +
  labs(x = "Parameters", y = "Relative importance (%)") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(expand = c(0, 0))

