# Load necessary libraries
library(tidyverse)
library(pROC)
library(ggplot2)
library(svglite)
library(mclust)

# Load pre-processed data
train_data <- read_tsv("Training_TAILVAR_score.txt")
test_data <- read_tsv("Testing_TAILVAR_score.txt")
stoploss_all <- read_tsv("stoploss_SNV_TAILVAR_score.txt")
stoploss_gnomAD <- stoploss_all %>% filter(gnomAD_exomes_POPMAX_AF > 0 | gnomAD_genomes_POPMAX_AF > 0)

# Compute ROC curve and AUROC for the train dataset
rf_roc_curve <- roc(train_data$Class, train_data$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the training set:", rf_auroc, "\n")

# Plot the distribution of TAILVAR Score (Train set: HGMD + ALFA)
plot1 <- ggplot(train_data, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Train set", color = "Train set") +
  scale_fill_manual(values = c("P" = "#EE9988", "B" = "#77AADD"), labels = c("P" = "HGMD DM", "B" = "ALFA AF>0.001")) +
  scale_color_manual(values = c("P" = "#BB4444", "B" = "#4477AA"), labels = c("P" = "HGMD DM", "B" = "ALFA AF>0.001")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "top", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
plot1
ggsave("Train_TAILVAR_distribtion.svg", plot = plot1, width = 6, height = 4)

# Compute ROC curve and AUROC for the test dataset
test_data$Class <- factor(test_data$Class, levels = c("VUS", "B", "P"))
test_data1 <- test_data %>% filter(Class %in% c("B","P"))
test_data2 <- test_data %>% filter(Class == "VUS")

rf_roc_curve <- roc(test_data1$Class, test_data1$TAILVAR, levels = c("B", "P"))
rf_auroc <- auc(rf_roc_curve)
cat("AUROC of the Random Forest model on the testing set:", rf_auroc, "\n")

# Plot the distribution of TAILVAR score (Test set: ClinVar + gnomAD)
plot2 <- ggplot(test_data1, aes(x = TAILVAR, fill = Class, color = Class)) +
  geom_histogram(binwidth = 0.05, size = 1.0, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "TAILVAR score", y = "Frequency", fill = "Test set", color = "Test set") +
  scale_fill_manual(values = c("P" = "#EE9988", "VUS" = "grey", "B" = "#77AADD"), labels = c("P" = "ClinVar P/LP", "VUS" = "VUS", "B" = "ClinVar B/LB + gnomAD AF>0.001")) +
  scale_color_manual(values = c("P" = "#BB4444", "VUS" = "darkgrey", "B" = "#4477AA"), labels = c("P" = "ClinVar P/LP", "VUS" = "VUS", "B" = "ClinVar B/LB + gnomAD AF>0.001")) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))+
  theme(legend.position = "top", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
plot2
ggsave("Test_TAILVAR_distribtion.svg", plot = plot2, width = 6, height = 4)

# Fit the Gaussian Mixture Model to all stoploss variants
gmm_model <- Mclust(stoploss_all$TAILVAR, G = 2)
gmm_means <- gmm_model$parameters$mean
gmm_sd <- sqrt(gmm_model$parameters$variance$sigmasq)
gmm_weights <- gmm_model$parameters$pro

pdf_1 <- function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * gmm_weights[1]
pdf_2 <- function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * gmm_weights[2]

likelihood_ratio_p <- function(x) pdf_2(x) / pdf_1(x)
x_vals <- seq(min(stoploss_all$TAILVAR), max(stoploss_all$TAILVAR), length.out = 1000)
lr_vals <- sapply(x_vals, likelihood_ratio_p)

jpeg("LR_TAILVAR_threshold.jpg", width = 6 * 300, height = 6 * 300, res = 300)
LR_plot <- plot(x_vals, lr_vals, type = "l", col = "blue", lwd = 2,
                xlab = "TAILVAR Score", ylab = "Likelihood Ratio",
                main = "Likelihood Ratio vs. TAILVAR Score")
abline(h = 9, col = "red", lty = 2)
abline(h = 1/9, col = "green", lty = 2)
dev.off()

# Find the threshold where LR = 9
threshold_p <- uniroot(function(x) likelihood_ratio_p(x) - 9, 
                       interval = c(0.80, 0.95))$root
threshold_b <- uniroot(function(x) likelihood_ratio_p(x) - 0.1111, 
                       interval = c(0.00, 0.90))$root
# Print the threshold
cat("Threshold TAILVAR score for LR = 9:", threshold_p, "\n")
cat("Threshold TAILVAR score for LR = 0.1111:", threshold_b, "\n")


# Plot the histogram with Gaussian curves and threshold
max_freq <- nrow(stoploss_all)/(21*2)
plot3 <- ggplot(stoploss_all, aes(x = TAILVAR)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "white", alpha = 0.7) +
  labs(x = "TAILVAR score", y = "Frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * max_freq, color = "#00B050", linetype = "solid", size = 1.5) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * max_freq, color = "#FEA04C", linetype = "solid", size = 1.5) +
  stat_function(fun = function(x) (dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) + dnorm(x, mean = gmm_means[2], sd = gmm_sd[2])) * max_freq, color = "black", linetype = "dashed", size = 1.0)
plot3
ggsave("stoploss_all_GMM.svg", plot = plot3, width = 6, height = 5)


# Fit the Gaussian Mixture Model to ClinVar VUS variants
gmm_model <- Mclust(test_data2$TAILVAR, G = 2)
gmm_means <- gmm_model$parameters$mean
gmm_sd <- sqrt(gmm_model$parameters$variance$sigmasq)
gmm_weights <- gmm_model$parameters$pro

pdf_1 <- function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * gmm_weights[1]
pdf_2 <- function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * gmm_weights[2]

max_freq <- nrow(test_data2)/(21*2)
plot4 <- ggplot(test_data2, aes(x = TAILVAR)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black", alpha = 0.7) +
  labs(x = "TAILVAR score", y = "Frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * max_freq,
                color = "#9DC3E6", linetype = "solid", size = 1.5) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * max_freq,
                color = "#FF75AD", linetype = "solid", size = 1.5) +
  geom_vline(xintercept = threshold_b, color = "#0070C0", linetype = "dashed", size = 1.0) +
  geom_vline(xintercept = threshold_p, color = "#C00000", linetype = "dashed", size = 1.0)
plot4
# Save the plot
ggsave("Clinvar_VUS_GMM.svg", plot = plot4, width = 6, height = 5)

# Plot the distribution of TAILVAR score (gnomAD_all) according to AF
gnomAD_AF <- c("gnomAD_exomes_AF", "gnomAD_exomes_POPMAX_AF", "gnomAD_genomes_AF", "gnomAD_genomes_POPMAX_AF")
stoploss_gnomAD$gnomAD_AF <- apply(stoploss_gnomAD[gnomAD_AF], 1, max, na.rm = TRUE)
stoploss_gnomAD <- stoploss_gnomAD %>% mutate(pAF = -log(gnomAD_AF))
plot5 <- ggplot(stoploss_gnomAD, aes(y = TAILVAR, x = pAF)) +
  geom_jitter(width = 0.5, height = 0.5, size = 1.0, alpha = 0.5, color = "grey") +
  geom_smooth(method = "auto", color = "red", size = 1.2, se = TRUE) +
  labs(title = "gnomAD stop-lost variants", y = "TAILVAR score", x = "-log (AF)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +  theme_classic() + 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 14, by = 2))
plot5
ggsave("gnomAD_AF_TAILVAR_relationship.svg", plot = plot5, width = 7, height = 4)

# Plot the relationship between extension length and TAILVAR score
stoploss_all <- stoploss_all %>% mutate(Extension_lengths = as.numeric(Extension_lengths)) %>% filter(is.na(Extension_lengths) == FALSE)
plot6 <- ggplot(stoploss_all, aes(y = TAILVAR, x = Extension_lengths)) +
  geom_jitter(width = 0.5, height = 0.5, size = 0.5, alpha = 0.1, color = "navy") +
  geom_smooth(method = "auto", color = "red", size = 1.2, se = TRUE) +
  labs(y = "TAILVAR score", x = "Extension_lengths") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +  theme_classic() + 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, by = 20))
plot6
ggsave("extension_length_TAILVAR_relationship.svg", plot = plot6, width = 6, height = 4)

# AUC-ROC plots
comp_scores <- c("CADD", "DANN", "FATHMM", "Eigen", "BayesDel_addAF", "BayesDel_noAF", "integrated_fitCons", "GERP", "phyloP100way", "phastCons100way")
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
  labs(title = "Train set: HGMD/ALFA", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values,")"))
auroc_plot1
ggsave("Train_dataset_AUROC.svg", plot = auroc_plot1, width = 6, height = 4)

auroc_plot2 <- ggroc(roc_list2, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "Test set: ClinVar/gnomAD", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values2,")"))
auroc_plot2
ggsave("Test_dataset_AUROC.svg", plot = auroc_plot2, width = 6, height = 4)

# Compare AUROC values
auroc_df_dev <- data.frame(Comparator = comparators, AUROC = auroc_values, Dataset = "HGMD/ALFA")
auroc_df_val <- data.frame(Comparator = comparators, AUROC = auroc_values2, Dataset = "ClinVar/gnomAD")
auroc_combined_df <- rbind(auroc_df_dev, auroc_df_val)
auroc_combined_df$Dataset <- factor(auroc_combined_df$Dataset, levels = c("HGMD/ALFA", "ClinVar/gnomAD"))

# Order the data frame by AUROC in decreasing order for the development dataset
auroc_combined_df$Comparator <- factor(auroc_combined_df$Comparator, levels = auroc_df_dev$Comparator[order(-auroc_df_dev$AUROC)])

# Replot with bars filled by Dataset
auroc_plot3 <- ggplot(auroc_combined_df, aes(x = Comparator, y = AUROC, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") + labs(x = "", y = "AUROC") +
  theme_classic() + theme(axis.title.y = element_text(size = 14), axis.text.x = element_text(angle = -45, hjust = 0), legend.position = "top") +
  scale_fill_manual(values = c("HGMD/ALFA" = "#C0A5DB", "ClinVar/gnomAD" = "#FFDF7F")) + coord_cartesian(ylim =c(0.5, 1))

auroc_plot3
ggsave("TAILVAR_AUROC_comparisons.svg", plot = auroc_plot3, width = 7, height = 4)
