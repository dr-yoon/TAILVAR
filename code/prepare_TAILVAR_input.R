# Load necessary libraries
library(tidyverse)

# Load preprocessed transcripts data
transcript_info <- read_tsv("mane_transcripts_info.txt")
transcript_info <- transcript_info %>% dplyr::select(ensembl_transcript_id, percentage_gene_gc_content, length_3utr, gc_3utr, stop_codons, next_stop_codon, ext_length, translate_tail)
transcript_info_colname <- c("ensembl_transcript_id", "percentage_gene_gc_content", "length_3utr", "gc_3utr", "ext_length", "translate_tail")
transcript_info_rename <- c("Feature", "Gene_GC", "UTR3_length", "UTR3_GC", "TailAA_counts", "TailAA_seq")

transcript_info <- transcript_info %>% rename_at(vars(transcript_info_colname), ~ transcript_info_rename)

# Load all possible stoplost SNV data
stoplost_all <- read_tsv("stoplost_SNV_filtered.txt")
stoplost_all <- stoplost_all %>% mutate(Class = "")
stoplost_all <- left_join(stoplost_all, transcript_info, by = c("Feature")) %>% distinct()
stoplost_all <- stoplost_all %>% mutate(TailAA_seq = paste0(Amino_acids, TailAA_seq))
stoplost_all <- stoplost_all %>% filter(is.na(TailAA_counts) == FALSE)

# filtered by gnomAD AF
stoplost_gnomAD <- stoplost_all %>% filter(gnomAD_exomes_POPMAX_AF > 0 | gnomAD_genomes_POPMAX_AF > 0)
gnomAD_MAF_0.001 <- stoplost_gnomAD %>% filter(gnomAD_exomes_POPMAX_AF > 0.001 | gnomAD_genomes_POPMAX_AF > 0.001)
gnomAD_MAF_0.001 <- gnomAD_MAF_0.001 %>% mutate(Class = "B")

# Load HGMD data
stoplost_HGMD <- read_tsv("stoplost_HGMD_filtered.txt")
stoplost_HGMD <- stoplost_HGMD %>% mutate(Class = "P")
stoplost_HGMD <- left_join(stoplost_HGMD, transcript_info, by = c("Feature")) %>% distinct()
stoplost_HGMD <- stoplost_HGMD %>% mutate(TailAA_seq = paste0(Amino_acids, TailAA_seq))
stoplost_HGMD <- stoplost_HGMD %>% filter(is.na(TailAA_counts) == FALSE)

# Training set: HGMD + gnomAD data
train_data <- rbind(stoplost_HGMD, gnomAD_MAF_0.001)

# filtered by EAS AF
ALFA_AF <- c("ALFA_European_AF", "ALFA_African_AF", "ALFA_Asian_AF", "ALFA_Other_AF", "ALFA_Total_AF")
ALFA_MAF_0.001 <- stoplost_all %>% filter(if_any(all_of(ALFA_AF), ~ . > 0.001))
ALFA_MAF_0.001 <- ALFA_MAF_0.001 %>% mutate(Class = "B")

# Load Clinvar data
stoplost_clinvar <- read_tsv("stoplost_Clinvar_filtered.txt")
stoplost_clinvar <- stoplost_clinvar %>% filter(Consequence %in% c("stop_lost")) %>%
  filter(ClinVar_CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic","Uncertain_significance","Benign", "Likely_benign", "Benign/Likely_benign", "Conflicting_classifications_of_pathogenicity")) %>%
  mutate(
    Class = case_when(str_detect(ClinVar_CLNSIG, regex("athogenic", ignore_case = TRUE)) ~ "P",
                               str_detect(ClinVar_CLNSIG, regex("ncertain", ignore_case = TRUE)) ~ "VUS", 
                               str_detect(ClinVar_CLNSIG, regex("enign", ignore_case = TRUE)) ~ "B",
                               str_detect(ClinVar_CLNSIG, regex("Conflicting", ignore_case = TRUE)) ~ "C",
                               TRUE ~ ClinVar_CLNSIG),
    )
stoplost_clinvar <- left_join(stoplost_clinvar, transcript_info, by = c("Feature")) %>% distinct()
stoplost_clinvar <- stoplost_clinvar %>% mutate(TailAA_seq = paste0(Amino_acids, TailAA_seq))

# Testing set: ClinVar + ALFA data
test_data <- rbind(stoplost_clinvar, ALFA_MAF_0.001)

# Add columns for each amino acid count
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (aa in amino_acid_columns) {train_data <- train_data %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {test_data <- test_data %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {stoplost_all <- stoplost_all %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {stoplost_gnomAD <- stoplost_gnomAD %>% mutate(!!aa := str_count(TailAA_seq, aa))}
hydrophobic_aa <- c("A", "I", "L", "M", "F", "V", "P")

train_data <- train_data %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/TailAA_counts)
test_data <- test_data %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/TailAA_counts)
stoplost_all <- stoplost_all %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/TailAA_counts)
stoplost_gnomAD <- stoplost_gnomAD %>% mutate(Hydrophobicity = rowSums(select(., all_of(hydrophobic_aa)))/TailAA_counts)

# Define features and target
comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_length", "UTR3_GC", "TailAA_counts")

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
  if (is.numeric(stoplost_all[[scores]])) {
    stoplost_all[[scores]][is.na(stoplost_all[[scores]])] <- median(stoplost_all[[scores]], na.rm = TRUE)
  }
}
for (scores in comp_scores) {
  if (is.numeric(stoplost_gnomAD[[scores]])) {
    stoplost_gnomAD[[scores]][is.na(stoplost_gnomAD[[scores]])] <- median(stoplost_gnomAD[[scores]], na.rm = TRUE)
  }
}

# Prepare the dataset
train_data <- train_data %>% filter(Class %in% c("B", "P"), is.na(TailAA_counts) == FALSE)
test_data <- test_data %>% filter(Class %in% c("B", "VUS", "P"), is.na(TailAA_counts) == FALSE)

# Write the filtered dataset to a file
write.table(train_data, "HGMD_gnomAD_train_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(test_data, "ClinVar_ALFA_test_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_gnomAD, "stoplost_gnomAD_prediction_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_all, "stoplost_SNV_prediction_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
