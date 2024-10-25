# Load necessary libraries
library(tidyverse)
library(pROC)
library(ggplot2)
library(svglite)
library(mclust)

# Load pre-processed data
train_data <- read_tsv("Training_TAILVAR_score.txt")
test_data <- read_tsv("Testing_TAILVAR_score.txt")
stoplost_all <- read_tsv("stoplost_SNV_TAILVAR_score.txt")
stoplost_gnomAD <- stoplost_all %>% filter(gnomAD_exomes_POPMAX_AF > 0 | gnomAD_genomes_POPMAX_AF > 0)

# Compute ROC curve and AUROC for the train dataset
rf_roc_curve <- roc(train_data$Class, train_data$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the training set:", rf_auroc, "\n")

# Plot the distribution of TAILVAR Score (HGMD + gnomAD)
plot1 <- ggplot(train_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P" = "#EE9988", "B" = "#77AADD"), labels = c("P" = "HGMD", "B" = "gnomAD")) +
  scale_color_manual(values = c("P" = "#BB4444", "B" = "#4477AA"), labels = c("P" = "HGMD", "B" = "gnomAD")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")
plot1
ggsave("Train_dataset_TAILVAR_distribtion.svg", plot = plot1, width = 6, height = 3)

# Compute ROC curve and AUROC for the test dataset
test_data$Class <- factor(test_data$Class, levels = c("VUS", "P", "B"))
test_data1 <- test_data %>% filter(Class %in% c("B","P"))
test_data2 <- test_data %>% filter(Class == "VUS")

rf_roc_curve <- roc(test_data1$Class, test_data1$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the testing set:", rf_auroc, "\n")

# Plot the distribution of TAILVAR score (ClinVar)
plot2 <- ggplot(test_data1, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Dataset", color = "Dataset") +
  scale_fill_manual(values = c("P" = "#EE9988", "VUS" = "grey", "B" = "#77AADD"), labels = c("P" = "P/LP", "VUS" = "VUS", "B" = "B/LB")) +
  scale_color_manual(values = c("P" = "#BB4444", "VUS" = "darkgrey", "B" = "#4477AA"), labels = c("P" = "P/LP", "VUS" = "VUS", "B" = "B/LB")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))+
  theme(legend.position = "none")
plot2
ggsave("Test_dataset_TAILVAR_distribtion.svg", plot = plot2, width = 6, height = 3)

# Create the plot with Gaussian curves scaled to match the histogram
gmm_model <- Mclust(test_data2$TAILVAR, G = 2)
gmm_means <- gmm_model$parameters$mean
gmm_sd <- sqrt(gmm_model$parameters$variance$sigmasq)
max_freq <- 20

plot3 <- ggplot(test_data2, aes(x = TAILVAR)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black", alpha = 0.7) +
  labs(x = "TAILVAR score", y = "Frequency") +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) +scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none") +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * max_freq,
                color = "#4477AA", linetype = "dashed", size = 1.0) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * max_freq,
                color = "#BB4444", linetype = "dashed", size = 1.0)
plot3 
ggsave("Clinvar_VUS_Gaussian_mixture.svg", plot = plot3, width = 6, height = 5)

# Plot the distribution of TAILVAR score (gnomAD_all) according to AF
gnomAD_AF <- c("gnomAD_exomes_AF", "gnomAD_exomes_POPMAX_AF", "gnomAD_genomes_AF", "gnomAD_genomes_POPMAX_AF")
stoplost_gnomAD$gnomAD_AF <- apply(stoplost_gnomAD[gnomAD_AF], 1, max, na.rm = TRUE)
stoplost_gnomAD <- stoplost_gnomAD %>% mutate(pAF = -log(gnomAD_AF))
plot4 <- ggplot(stoplost_gnomAD, aes(y = TAILVAR, x = pAF)) +
  geom_point(size = 1, alpha = 0.9, color = "grey") +
  geom_smooth(method = "auto", color = "red", size = 1.2, se = TRUE) +
  labs(title = "gnomAD allele frequency", y = "TAILVAR score", x = "-log (AF)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
plot4
ggsave("gnomAD_AF_TAILVAR_relationship.svg", plot = plot4, width = 6, height = 5)

# AUC-ROC plots
comparators <- c(all_of(comp_scores),"TAILVAR")
AUROC_data <- train_data
AUROC_data$Class <- ifelse(AUROC_data$Class == "P", 1, 0)

AUROC_data2 <- test_data1
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
auroc_plot1 <- ggroc(roc_list, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "ROC curves (HGMD/gnomAD)", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values,")"))
auroc_plot1
ggsave("Train_dataset_TAILVAR_AUROC.svg", plot = auroc_plot1, width = 8, height = 6)

auroc_plot2 <- ggroc(roc_list2, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "ROC curves (ClinVar)", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right") + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values2,")"))
auroc_plot2
ggsave("Validation_dataset_TAILVAR_AUROC.svg", plot = auroc_plot2, width = 8, height = 6)

# Compare AUROC values
auroc_df_dev <- data.frame(Comparator = comparators, AUROC = auroc_values, Dataset = "HGMD/gnomAD")
auroc_df_val <- data.frame(Comparator = comparators, AUROC = auroc_values2, Dataset = "ClinVar")
auroc_combined_df <- rbind(auroc_df_dev, auroc_df_val)
auroc_combined_df$Dataset <- factor(auroc_combined_df$Dataset, levels = c("HGMD/gnomAD", "ClinVar"))

# Order the data frame by AUROC in decreasing order for the development dataset
auroc_combined_df$Comparator <- factor(auroc_combined_df$Comparator, levels = auroc_df_dev$Comparator[order(-auroc_df_dev$AUROC)])

# Replot with bars filled by Dataset
auroc_plot3 <- ggplot(auroc_combined_df, aes(x = Comparator, y = AUROC, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") + labs(x = "", y = "AUROC") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("HGMD/gnomAD" = "skyblue", "ClinVar" = "orange")) + coord_cartesian(ylim =c(0.5, 1))
auroc_plot3
ggsave("TAILVAR_AUROC_comparisons.svg", plot = auroc_plot3, width = 10, height = 6)
