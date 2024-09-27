# Load necessary libraries
library(tidyverse)
library(pROC)

# Load pre-processed data
model_data <- read_tsv("Development_dataset_TAILVAR_score.txt")
validation_data <- read_tsv("Validation_dataset_TAILVAR_score.txt")
stoplost_all <- read_tsv("stoplost_SNV_TAILVAR_score.txt")
stoplost_gnomAD <- read_tsv("stoplost_gnomAD_TAILVAR_score.txt")

# Define feature groups for model input
comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", 
                 "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_GC", "UTR3_length", "TailAA_counts")
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Compute ROC curve and AUROC for the HGMD+gnomAD dataset
rf_roc_curve <- roc(model_data$Class, model_data$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the test dataset:", rf_auroc, "\n")

# Plot the distribution of TAILVAR Score (HGMD + gnomAD)
ggplot(model_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P" = "#EE9988", "B" = "#77AADD"), labels = c("P" = "HGMD", "B" = "gnomAD")) +
  scale_color_manual(values = c("P" = "#BB4444", "B" = "#4477AA"), labels = c("P" = "HGMD", "B" = "gnomAD")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")

# Compute ROC curve and AUROC for the ClinVar dataset
rf_roc_curve <- roc(validation_data$Class, validation_data$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the validation dataset:", rf_auroc, "\n")

# Plot the distribution of TAILVAR score (ClinVar)
validation_data$Class <- factor(validation_data$Class, levels = c("VUS", "P", "B"))
ggplot(validation_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P" = "#EE9988", "VUS" = "grey", "B" = "#77AADD"), labels = c("P" = "P/LP", "VUS" = "VUS", "B" = "B/LB")) +
  scale_color_manual(values = c("P" = "#BB4444", "VUS" = "darkgrey", "B" = "#4477AA"), labels = c("P" = "P/LP", "VUS" = "VUS", "B" = "B/LB")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))+
  theme(legend.position = "none")

# Plot the distribution of TAILVAR score (gnomAD_all) according to AF
stoplost_gnomAD <- stoplost_gnomAD %>% mutate(pAF = -log(gnomAD_AF))
ggplot(stoplost_gnomAD, aes(y = TAILVAR, x = pAF)) +
  geom_point(size = 1, alpha = 0.9, color = "grey") +
  geom_smooth(method = "auto", color = "red", size = 1.2, se = TRUE) +
  labs(title = "gnomAD allele frequency", y = "TAILVAR score", x = "-log (AF)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

# AUC-ROC plot
# Define the features and the labels
comparators <- c(all_of(comp_scores),"TAILVAR")
AUROC_data <- model_data
AUROC_data$Class <- ifelse(AUROC_data$Class == "P", 1, 0)

AUROC_data2 <- validation_data %>% filter(Class %in% c("B","P"))
AUROC_data2$Class <- ifelse(AUROC_data2$Class == "P", 1, 0)

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
  labs(title = "ROC curves (HGMD/gnomAD)", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values,")"))

ggroc(roc_list2, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "ROC curves (ClinVar)", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values2,")"))

# Compare AUROC values
auroc_df_dev <- data.frame(Comparator = comparators, AUROC = auroc_values, Dataset = "HGMD/gnomAD")
auroc_df_val <- data.frame(Comparator = comparators, AUROC = auroc_values2, Dataset = "ClinVar")
auroc_combined_df <- rbind(auroc_df_dev, auroc_df_val)
auroc_combined_df$Dataset <- factor(auroc_combined_df$Dataset, levels = c("HGMD/gnomAD", "ClinVar"))

# Order the data frame by AUROC in decreasing order for the development dataset
auroc_combined_df$Comparator <- factor(auroc_combined_df$Comparator, levels = auroc_df_dev$Comparator[order(-auroc_df_dev$AUROC)])

# Replot with bars filled by Dataset
ggplot(auroc_combined_df, aes(x = Comparator, y = AUROC, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") + labs(x = "", y = "AUROC") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("HGMD/gnomAD" = "skyblue", "ClinVar" = "orange")) + coord_cartesian(ylim =c(0.5, 1))

