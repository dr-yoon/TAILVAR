# Load necessary libraries
library(tidyverse)
library(caret)
library(randomForest)


# Load pre-processed data
model_data <- read_tsv("HGMD_gnomAD_model_data.txt")
validation_data <- read_tsv("ClinVar_validation_data.txt")
stoplost_all <- read_tsv("stoplost_SNV_prediction_data.txt")
stoplost_gnomAD <- read_tsv("stoplost_gnomAD_prediction_data.txt")


comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_GC", "UTR3_length", "TailAA_counts")
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# RF model development
model_input <- model_data %>% dplyr::select("Class", all_of(comp_scores), all_of(gene_feature), all_of(amino_acid_columns)) %>%
  mutate(Class = factor(Class, levels = c("B", "P")))
set.seed(123)
trainIndex <- createDataPartition(model_input$Class, p = 0.8, list = FALSE)
train_data <- model_input[trainIndex,]
test_data <- model_input[-trainIndex,]

# Set up cross-validation
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Optimize hyper-parameters: mtry, ntree
tune_grid <- expand.grid(mtry = c(2, 5, 10, 20, 30))
ntree_values <- c(10, 25, 50, 100, 150, 200)
results <- data.frame()

for (ntree in ntree_values) {
  rf_tune <- train(Class ~ ., data = train_data, method = "rf", 
                   trControl = train_control, tuneGrid = tune_grid, metric = "ROC", 
                   ntree = ntree)
  for (mtry in tune_grid$mtry) {
    roc_value <- rf_tune$results$ROC[rf_tune$results$mtry == mtry]
    results <- rbind(results, data.frame(ntree = ntree, mtry = mtry, ROC = roc_value))
  }
}

# Plot the AUROC vs. ntree for each mtry value
ggplot(results, aes(x = ntree, y = ROC, color = factor(mtry))) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "AUROC vs. hyper-parameters (ntree, mtry)",
       x = "Number of Trees (ntree)", y = "AUROC", color = "mtry") +
  theme_minimal() + theme(legend.position = "right")

# Select the best combination of mtry and ntree
selected_params <- c(10, 100)
cat("Selected mtry:", selected_params[1], "with ntree:", selected_params[2], "\n")

# Train the final model with the best parameters
final_rf_model <- randomForest(Class ~ ., data = train_data, 
                               mtry = selected_params[1], ntree = selected_params[2])

# Predict probabilities for the overall dataset and name the score as TAILVAR_score
model_data$TAILVAR <- predict(final_rf_model, model_data, type = "prob")[, "P"]
validation_data$TAILVAR <- predict(final_rf_model, validation_data, type = "prob")[, "P"]
stoplost_all$TAILVAR <- predict(final_rf_model, stoplost_all, type = "prob")[, "P"]
stoplost_gnomAD$TAILVAR <- predict(final_rf_model, stoplost_gnomAD, type = "prob")[, "P"]

write.table(model_data, "Development_dataset_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(validation_data, "Validation_dataset_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_gnomAD, "stoplost_gnomAD_TAILVAR_score.", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_all, "stoplost_SNV_TAILVAR_score.txt", row.names = FALSE, sep = "\t", quote = FALSE)