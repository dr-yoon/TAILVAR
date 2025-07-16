# Load necessary libraries
library(tidyverse)
library(pROC)
library(ggplot2)
library(mclust)

# Load pre-processed data
train_data <- read_tsv("Train_TAILVAR_score.txt")
test_data <- read_tsv("Test_TAILVAR_score.txt")
SNV_data <- read_tsv("stoploss_SNV_TAILVAR_score.txt")

# Compute ROC curve and AUROC for the train dataset
rf_roc_curve <- roc(train_data$Class, train_data$TAILVAR, levels = c("B_LB", "P_LP"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the training set:", rf_auroc, "\n")

# Plot the distribution of TAILVAR Score (HGMD + gnomAD)
ggplot(train_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 0.8, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", title = "Training set", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P_LP" = "#EE9988", "B_LB" = "#77AADD"), labels = c("P_LP" = "P/LP", "B_LB" = "B/LB")) +
  scale_color_manual(values = c("P_LP" = "#BB4444", "B_LB" = "#4477AA"), labels = c("P_LP" = "P/LP", "B_LB" = "B/LB")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "right", legend.title = element_blank())

# Compute ROC curve and AUROC for the test dataset
rf_roc_curve2 <- roc(test_data$Class, test_data$TAILVAR, levels = c("B_LB", "P_LP"))
rf_auroc2 <- auc(rf_roc_curve2)
cat("AUROC of the Random Forest model on the test set:", rf_auroc2, "\n")

# Plot the distribution of TAILVAR Score (Clinvar/ALFA/RME/AOU)
ggplot(test_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 0.8, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", title = "Test set", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P_LP" = "#EE9988", "B_LB" = "#77AADD"), labels = c("P_LP" = "P/LP", "B_LB" = "B/LB")) +
  scale_color_manual(values = c("P_LP" = "#BB4444", "B_LB" = "#4477AA"), labels = c("P_LP" = "P/LP", "B_LB" = "B/LB")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "right", legend.title = element_blank())

# AUC-ROC plot
# Define the features and the labels
comp_scores <- c("CADD", "DANN", "Eigen_PC", "BayesDel_noAF", "BayesDel_addAF", "FATHMM_MKL", "fitCons_int", "MutationTaster", "VEST4", "GPN_MSA", "GERP", "phyloP100", "phastCons100")
comparators <- c(all_of(comp_scores),"TAILVAR")
AUROC_data <- train_data
AUROC_data$Class <- ifelse(AUROC_data$Class == "P_LP", 1, 0)

AUROC_data2 <- test_data %>% filter(Class %in% c("P_LP", "B_LB"))
AUROC_data2$Class <- ifelse(AUROC_data2$Class == "P_LP", 1, 0)

roc_list <- list()
auroc_values <- c()

roc_list2 <- list()
auroc_values2 <- c()

# Calculate ROC for each feature and store the AUROC values
for (score in comparators) {
  roc_obj <- roc(AUROC_data$Class, AUROC_data[[score]], ci=TRUE)
  roc_list[[score]] <- roc_obj
  auroc_values <- c(auroc_values, round(auc(roc_obj), 3))
}

for (score in comparators) {
  roc_obj2 <- roc(AUROC_data2$Class, AUROC_data2[[score]], ci=TRUE)
  roc_list2[[score]] <- roc_obj2
  auroc_values2 <- c(auroc_values2, round(auc(roc_obj2), 3))
}

# Plot AUROC comparisons
ggroc(roc_list, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "Training set", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values,")"))

ggroc(roc_list2, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "Test set", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values2,")"))

# Compare AUROC values
auroc_df_dev <- data.frame(Comparator = comparators, AUROC = auroc_values, Dataset = "Training")
auroc_df_val <- data.frame(Comparator = comparators, AUROC = auroc_values2, Dataset = "Test")
auroc_combined_df <- rbind(auroc_df_dev, auroc_df_val)
auroc_combined_df$Dataset <- factor(auroc_combined_df$Dataset, levels = c("Training", "Test"))

# Order the data frame by AUROC in decreasing order for the development dataset
auroc_combined_df$Comparator <- factor(auroc_combined_df$Comparator, levels = auroc_df_dev$Comparator[order(-auroc_df_dev$AUROC)])

# Replot with bars filled by Dataset
ggplot(auroc_combined_df, aes(x = Comparator, y = AUROC, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") + labs(x = "", y = "AUROC") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.title = element_blank()) +
  scale_fill_manual(values = c("Training" = "#C0A5DB", "Test" = "#FFDF7F")) + coord_cartesian(ylim =c(0.5, 1))
