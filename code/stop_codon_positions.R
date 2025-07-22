# Load necessary libraries
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]

MANE_transcripts <- read_tsv("MANE_stop_codon_info.tsv") %>% mutate(ensembl_transcript_id = sub("\\..*$", "", transcript_id)) %>% dplyr::select(-transcript_id)

# Identify splitted transcripts
MANE_transcripts_count <- MANE_transcripts %>% group_by(ensembl_transcript_id) %>% summarize(count = n(), .groups = "drop")
no_splitted_transcript <- MANE_transcripts_count %>% filter(count == 1)
splitted_transcript <- MANE_transcripts_count %>% filter(count == 2)

# Identify stop codon positions for non-splitted transcripts
no_splitted_transcript <- MANE_transcripts %>% filter(ensembl_transcript_id %in% no_splitted_transcript$ensembl_transcript_id) %>%
  mutate(start = as.integer(start), end = as.integer(end),
    POS_1 = if_else(strand == "+", start, end), POS_2 = if_else(strand == "+", start + 1, end - 1), POS_3 = if_else(strand == "+", start + 2, end - 2)) %>%
  dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Identify stop codon positions for splitted transcripts
splitted_transcript <- MANE_transcripts %>% filter(ensembl_transcript_id %in% splitted_transcript$ensembl_transcript_id) %>% group_by(ensembl_transcript_id) %>% arrange(start, .by_group = TRUE) %>% mutate(row_num = row_number()) %>%
  pivot_wider(id_cols = c(ensembl_transcript_id, strand), names_from = row_num, values_from = c(start, end), names_glue = "POS{row_num}_{.value}") %>%
  mutate(start = as.integer(start), end = as.integer(end),
    POS_1 = if_else(strand == "+", POS1_start, POS2_end),
    POS_2 = if_else(strand == "+", if_else(POS1_start == POS1_end, POS2_start, POS1_end), if_else(POS2_start == POS2_end, POS1_end, POS2_start)),
    POS_3 = if_else(strand == "+", POS2_end, POS1_start)) %>% dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Combine non-splitted and splitted transcript positions
stop_codon_pos <- bind_rows(no_splitted_transcript, splitted_transcript)

# Merge transcript information with stop codon positions
MANE_transcripts <- inner_join(MANE_transcripts, stop_codon_pos, by = "ensembl_transcript_id") %>% distinct(ensembl_transcript_id, .keep_all = TRUE) %>% dplyr::select(-start, -end)

# Extract variant position and convert to a VCF data format
vcf_data_list <- list(
  MANE_transcripts %>% dplyr::select(chromosome, POS_1),
  MANE_transcripts %>% dplyr::select(chromosome, POS_2),
  MANE_transcripts %>% dplyr::select(chromosome, POS_3)
)

vcf_data <- bind_rows(lapply(vcf_data_list, function(df) {
  names(df) <- c("chromosome", "POS")
  df %>% rowwise()})) %>% arrange(chromosome, POS) %>% unique()

write.table(vcf_data, paste0(file_name,".tsv"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
