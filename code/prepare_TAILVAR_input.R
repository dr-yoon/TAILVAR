# Load necessary libraries
library(tidyverse)

# Load preprocessed transcripts data
transcript_info <- read_tsv("mane_transcripts_info.txt")
transcript_info <- transcript_info %>% dplyr::select(ensembl_transcript_id, percentage_gene_gc_content, length_3utr, gc_3utr, stop_codons, next_stop_codon, ext_length, translate_tail)
transcript_info_colname <- c("ensembl_transcript_id", "percentage_gene_gc_content", "length_3utr", "gc_3utr", "ext_length", "translate_tail")
transcript_info_rename <- c("Feature", "Gene_GC", "UTR3_length", "UTR3_GC", "Extension_lengths", "Extension_seq")

transcript_info <- transcript_info %>% rename_at(vars(transcript_info_colname), ~ transcript_info_rename)

# Load all possible stoploss SNV data
stoploss_all <- read_tsv("stoploss_SNV_filtered.txt")
stoploss_all <- stoploss_all %>% mutate(Class = "")
stoploss_all <- left_join(stoploss_all, transcript_info, by = c("Feature")) %>% distinct()
stoploss_all <- stoploss_all %>% mutate(Extension_seq = paste0(Amino_acids, Extension_seq))
stoploss_all <- stoploss_all %>% filter(is.na(Extension_lengths) == FALSE)

# Load Clinvar, HGMD data
stoploss_clinvar <- read_tsv("stoploss_Clinvar_filtered.txt")
stoploss_clinvar <- stoploss_clinvar %>% filter(Consequence %in% c("stop_lost")) %>%
  filter(ClinVar_CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic","Uncertain_significance","Benign", "Likely_benign", "Benign/Likely_benign")) %>%
  mutate(
    Class = case_when(str_detect(ClinVar_CLNSIG, regex("athogenic", ignore_case = TRUE)) ~ "P",
                      str_detect(ClinVar_CLNSIG, regex("ncertain", ignore_case = TRUE)) ~ "VUS", 
                      str_detect(ClinVar_CLNSIG, regex("enign", ignore_case = TRUE)) ~ "B",
                      str_detect(ClinVar_CLNSIG, regex("Conflicting", ignore_case = TRUE)) ~ "C",
                      TRUE ~ ClinVar_CLNSIG),
  )
stoploss_clinvar <- left_join(stoploss_clinvar, transcript_info, by = c("Feature")) %>% distinct()
stoploss_clinvar <- stoploss_clinvar %>% mutate(Extension_seq = paste0(Amino_acids, Extension_seq)) %>% filter(is.na(Extension_lengths) == FALSE)

stoploss_HGMD <- read_tsv("stoploss_HGMD_filtered.txt")
stoploss_HGMD <- stoploss_HGMD %>% mutate(Class = "P")
stoploss_HGMD <- left_join(stoploss_HGMD, transcript_info, by = c("Feature")) %>% distinct()
stoploss_HGMD <- stoploss_HGMD %>% mutate(Extension_seq = paste0(Amino_acids, Extension_seq))
stoploss_HGMD <- stoploss_HGMD %>% filter(is.na(Extension_lengths) == FALSE)

# Load Population DB: gnomAD, ALFA > 0.001
gnomAD_MAF_0.001 <- stoploss_all %>% filter(gnomAD_exomes_POPMAX_AF > 0.001 | gnomAD_genomes_POPMAX_AF > 0.001)
gnomAD_MAF_0.001 <- gnomAD_MAF_0.001 %>% mutate(Class = "B")

ALFA_AF <- c("ALFA_European_AF", "ALFA_African_AF", "ALFA_Asian_AF", "ALFA_Other_AF", "ALFA_Total_AF")
ALFA_MAF_0.001 <- stoploss_all %>% filter(if_any(all_of(ALFA_AF), ~ . > 0.001))
ALFA_MAF_0.001 <- ALFA_MAF_0.001 %>% mutate(Class = "B")


# Training set: HGMD + ALFA data
train_data <- rbind(stoploss_HGMD, ALFA_MAF_0.001) %>% unique()

# Testing set: ClinVar + gnomAD data
test_data <- rbind(stoploss_clinvar, gnomAD_MAF_0.001) %>% unique()

# Add columns for each amino acid count
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (aa in amino_acid_columns) {train_data <- train_data %>% mutate(!!aa := str_count(Extension_seq, aa))}
for (aa in amino_acid_columns) {test_data <- test_data %>% mutate(!!aa := str_count(Extension_seq, aa))}
for (aa in amino_acid_columns) {stoploss_all <- stoploss_all %>% mutate(!!aa := str_count(Extension_seq, aa))}

hydrophobic_aa <- c("A", "C", "I", "L", "M", "F", "V")
train_data <- train_data %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/Extension_lengths)
test_data <- test_data %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/Extension_lengths)
stoploss_all <- stoploss_all %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/Extension_lengths)

# Define features and target
comp_scores <- c("CADD", "DANN", "FATHMM", "Eigen", "BayesDel_addAF", "BayesDel_noAF", "integrated_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_length", "UTR3_GC", "Extension_lengths")

# Impute missing values with the median (for numeric data)
for (scores in comp_scores) {
  if (is.numeric(train_data[[scores]])) {
    train_data[[scores]][is.na(train_data[[scores]])] <- median(train_data[[scores]], na.rm = TRUE)
  }
}
for (scores in comp_scores) {
  if (is.numeric(test_data[[scores]])) {
    test_data[[scores]][is.na(test_data[[scores]])] <- median(test_data[[scores]], na.rm = TRUE)
  }
}
for (scores in comp_scores) {
  if (is.numeric(stoploss_all[[scores]])) {
    stoploss_all[[scores]][is.na(stoploss_all[[scores]])] <- median(stoploss_all[[scores]], na.rm = TRUE)
  }
}

# Prepare the dataset
train_data <- train_data %>% filter(Class %in% c("B", "P"), is.na(Extension_lengths) == FALSE)
test_data <- test_data %>% filter(Class %in% c("B", "VUS", "P"), is.na(Extension_lengths) == FALSE)

# Write the filtered dataset to a file
write.table(train_data, "Train_data_input.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(test_data, "Test_data_input.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoploss_all, "stoploss_SNV_input.txt", row.names = FALSE, sep = "\t", quote = FALSE)
