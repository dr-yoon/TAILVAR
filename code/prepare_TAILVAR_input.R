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
Pop_AF <- c("gnomAD_exomes_AF", "gnomAD_exomes_POPMAX_AF", "gnomAD_genomes_AF", "gnomAD_genomes_POPMAX_AF")
stoplost_all$gnomAD_AF <- apply(stoplost_all[Pop_AF], 1, max, na.rm = TRUE)

stoplost_gnomAD <- stoplost_all %>% filter(gnomAD_AF > 0)
gnomAD_MAF_0.001 <- stoplost_gnomAD %>% filter(gnomAD_AF > 0.001)
gnomAD_MAF_0.001 <- gnomAD_MAF_0.001 %>% mutate(Class = "B")

# Load HGMD stoplost data
stoplost_HGMD <- read_tsv("stoplost_HGMD_filtered.txt")
stoplost_HGMD <- stoplost_HGMD %>% mutate(Class = "P")
stoplost_HGMD <- left_join(stoplost_HGMD, transcript_info, by = c("Feature")) %>% distinct()
stoplost_HGMD <- stoplost_HGMD %>% mutate(TailAA_seq = paste0(Amino_acids, TailAA_seq))
stoplost_HGMD <- stoplost_HGMD %>% filter(is.na(TailAA_counts) == FALSE)
stoplost_HGMD$gnomAD_AF <- apply(stoplost_HGMD[Pop_AF], 1, max, na.rm = TRUE)

# combined HGMD and gnomAD stoplost data
model_data <- rbind(stoplost_HGMD, gnomAD_MAF_0.001)

validation_data <- read_tsv("stoplost_Clinvar_filtered.txt")
validation_data <- validation_data %>% filter(Consequence %in% c("stop_lost")) %>%
  filter(ClinVar_CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic","Uncertain_significance","Benign", "Likely_benign", "Benign/Likely_benign", "Conflicting_classifications_of_pathogenicity")) %>%
  mutate(
    Class = case_when(str_detect(ClinVar_CLNSIG, regex("athogenic", ignore_case = TRUE)) ~ "P",
                               str_detect(ClinVar_CLNSIG, regex("ncertain", ignore_case = TRUE)) ~ "VUS", 
                               str_detect(ClinVar_CLNSIG, regex("enign", ignore_case = TRUE)) ~ "B",
                               str_detect(ClinVar_CLNSIG, regex("Conflicting", ignore_case = TRUE)) ~ "C",
                               TRUE ~ ClinVar_CLNSIG),
    )
validation_data <- left_join(validation_data, transcript_info, by = c("Feature")) %>% distinct()
validation_data <- validation_data %>% mutate(TailAA_seq = paste0(Amino_acids, TailAA_seq))
validation_data$gnomAD_AF <- apply(validation_data[Pop_AF], 1, max, na.rm = TRUE)

# Add columns for each amino acid count
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (aa in amino_acid_columns) {model_data <- model_data %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {validation_data <- validation_data %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {stoplost_all <- stoplost_all %>% mutate(!!aa := str_count(TailAA_seq, aa))}
for (aa in amino_acid_columns) {stoplost_gnomAD <- stoplost_gnomAD %>% mutate(!!aa := str_count(TailAA_seq, aa))}

# Define features and target
comp_scores <- c("CADD", "DANN", "FATHMM", "EIGEN", "BayesDel_addAF", "BayesDel_noAF", "int_fitCons", "GERP", "phyloP100way", "phastCons100way")
gene_feature <- c("Gene_GC", "UTR3_length", "UTR3_GC", "TailAA_counts")

# Impute missing values with the median (for numeric data)
for (scores in comp_scores) {
  if (is.numeric(model_data[[scores]])) {
    model_data[[scores]][is.na(model_data[[scores]])] <- median(model_data[[scores]], na.rm = TRUE)
  }
}
for (scores in comp_scores) {
  if (is.numeric(validation_data[[scores]])) {
    validation_data[[scores]][is.na(validation_data[[scores]])] <- median(validation_data[[scores]], na.rm = TRUE)
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
model_data <- model_data %>% filter(Class %in% c("B", "P"), is.na(TailAA_counts) == FALSE)
validation_data <- validation_data %>% filter(Class %in% c("B", "VUS", "P"), is.na(TailAA_counts) == FALSE)

# Write the filtered dataset to a file
write.table(model_data, "HGMD_gnomAD_model_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(validation_data, "ClinVar_validation_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_gnomAD, "stoplost_gnomAD_prediction_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(stoplost_all, "stoplost_SNV_prediction_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)