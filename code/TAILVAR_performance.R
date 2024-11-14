# Load pre-processed data
train_data <- read_tsv("Training_TAILVAR_score.txt")
test_data <- read_tsv("Testing_TAILVAR_score.txt")
stoplost_all <- read_tsv("Stoplost_SNV_TAILVAR_score.txt")
stoplost_gnomAD <- stoplost_all %>% filter(gnomAD_exomes_POPMAX_AF > 0 | gnomAD_genomes_POPMAX_AF > 0)

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
ggsave("Test_TAILVAR_distribtion.svg", plot = plot2, width = 6, height = 4)

# Fit the Gaussian Mixture Model
gmm_model <- Mclust(test_data2$TAILVAR, G = 2)
gmm_means <- gmm_model$parameters$mean
gmm_sd <- sqrt(gmm_model$parameters$variance$sigmasq)
gmm_weights <- gmm_model$parameters$pro

# Define the PDFs for each component
pdf_1 <- function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * gmm_weights[1]
pdf_2 <- function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * gmm_weights[2]

# Define the likelihood ratio function
likelihood_ratio_p <- function(x) pdf_2(x) / pdf_1(x) # Assuming component 2 is P/LP
x_vals <- seq(min(test_data2$TAILVAR), max(test_data2$TAILVAR), length.out = 1000)
lr_vals <- sapply(x_vals, likelihood_ratio)
# Find the threshold where LR = 9
threshold_p <- uniroot(function(x) likelihood_ratio_p(x) - 9, 
                     interval = c(0.80, 0.99))$root
cat("Threshold TAILVAR score for LR = 9:", threshold_p, "\n")

# Plot the histogram with Gaussian curves and threshold
max_freq <- 15
plot3 <- ggplot(test_data2, aes(x = TAILVAR)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black", alpha = 0.7) +
  labs(x = "TAILVAR score", y = "Frequency") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[1], sd = gmm_sd[1]) * max_freq,
                color = "#4477AA", linetype = "dashed", size = 1.0) +
  stat_function(fun = function(x) dnorm(x, mean = gmm_means[2], sd = gmm_sd[2]) * max_freq,
                color = "#BB4444", linetype = "dashed", size = 1.0) +
  geom_vline(xintercept = threshold, color = "red", linetype = "solid", size = 1.0) +
  annotate("text", x = threshold_p, y = max_freq * 0.9, label = paste0("Threshold: ", round(threshold, 2)), 
           color = "red", size = 5, hjust = -0.1)
ggsave("Clinvar_VUS_Gaussian_mixture_with_threshold.svg", plot = plot3, width = 6, height = 5)

# Plot the distribution of TAILVAR score (gnomAD_all) according to AF
gnomAD_AF <- c("gnomAD_exomes_AF", "gnomAD_exomes_POPMAX_AF", "gnomAD_genomes_AF", "gnomAD_genomes_POPMAX_AF")
stoplost_gnomAD$gnomAD_AF <- apply(stoplost_gnomAD[gnomAD_AF], 1, max, na.rm = TRUE)
stoplost_gnomAD <- stoplost_gnomAD %>% mutate(pAF = -log(gnomAD_AF))
plot4 <- ggplot(stoplost_gnomAD, aes(y = TAILVAR, x = pAF)) +
  geom_jitter(width = 0.5, height = 0.5, size = 1.0, alpha = 0.5, color = "grey") +
  geom_smooth(method = "auto", color = "red", size = 1.2, se = TRUE) +
  labs(title = "gnomAD stop-lost variants", y = "TAILVAR score", x = "-log (AF)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +  theme_classic() + 
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 14, by = 2))
ggsave("gnomAD_AF_TAILVAR_relationship.svg", plot = plot4, width = 6, height = 4)

# AUC-ROC plots
comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
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
ggsave("Train_dataset_AUROC.svg", plot = auroc_plot1, width = 6, height = 4)

auroc_plot2 <- ggroc(roc_list2, legacy.axes = TRUE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = "Test set: ClinVar/gnomAD", x = "1 - Specificity", y = "Sensitivity") +
  theme_classic() + theme(legend.position = "right", axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + scale_color_discrete(name = "Comparators", labels = paste(comparators, "(AUROC =", auroc_values2,")"))
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
  theme_classic() + theme(axis.title.y = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("HGMD/ALFA" = "#C0A5DB", "ClinVar/gnomAD" = "#FFDF7F")) + coord_cartesian(ylim =c(0.5, 1))
ggsave("TAILVAR_AUROC_comparisons.svg", plot = auroc_plot3, width = 6, height = 4)
