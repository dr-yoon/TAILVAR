# Load necessary libraries
library(tidyverse)
library(biomaRt)
library(stringr)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)
vcf_body=args[1]

MANE_transcripts <- read_tsv("MANE_stop_codon_info.tsv") %>% mutate(ensembl_transcript_id = sub("\\..*$", "", transcript_id)) %>% dplyr::select(-transcript_id)

# Obtain coding and 3'UTR data from the Ensembl database
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#transcript_coding <- getBM(attributes = c("ensembl_transcript_id", "coding"), filters = "ensembl_transcript_id", values = MANE_transcripts$ensembl_transcript_id, mart = ensembl)
#transcript_3utr <- getBM(attributes = c("ensembl_transcript_id", "3_utr_start", "3_utr_end", "strand", "3utr"), filters = "ensembl_transcript_id", values = MANE_transcripts$ensembl_transcript_id, mart = ensembl)

# Read pre-downloaded data: since querying the Ensembl database can be time-consuming, we use pre-downloaded data files
transcript_coding <- read_tsv("MANE_transcripts_coding.txt")
transcript_3utr <- read_tsv("MANE_transcripts_3utr.txt")

# Identify splitted transcripts
MANE_transcripts_count <- MANE_transcripts %>% group_by(ensembl_transcript_id) %>% summarize(count = n(), .groups = "drop")
no_splitted_transcript <- MANE_transcripts_count %>% filter(count == 1)
splitted_transcript <- MANE_transcripts_count %>% filter(count == 2)

# Identify stop codon positions for non-splitted transcripts
no_splitted_transcript <- MANE_transcripts %>% filter(ensembl_transcript_id %in% no_splitted_transcript$ensembl_transcript_id) %>%
  mutate(POS_1 = if_else(strand == "+", start, end), POS_2 = if_else(strand == "+", start + 1, end - 1), POS_3 = if_else(strand == "+", start + 2, end - 2)) %>%
  dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Identify stop codon positions for splitted transcripts
splitted_transcript <- MANE_transcripts %>% filter(ensembl_transcript_id %in% splitted_transcript$ensembl_transcript_id) %>% group_by(ensembl_transcript_id) %>% arrange(start, .by_group = TRUE) %>% mutate(row_num = row_number()) %>%
  pivot_wider(id_cols = c(ensembl_transcript_id, strand), names_from = row_num, values_from = c(start, end), names_glue = "POS{row_num}_{.value}") %>%
  mutate(
    POS_1 = if_else(strand == "+", POS1_start, POS2_end),
    POS_2 = if_else(strand == "+", if_else(POS1_start == POS1_end, POS2_start, POS1_end), if_else(POS2_start == POS2_end, POS1_end, POS2_start)),
    POS_3 = if_else(strand == "+", POS2_end, POS1_start)) %>% dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Combine non-splitted and splitted transcript positions
stop_codon_pos <- bind_rows(no_splitted_transcript, splitted_transcript)

# Merge transcript information with stop codon positions
MANE_transcripts <- inner_join(MANE_transcripts, stop_codon_pos, by = "ensembl_transcript_id") %>% distinct(ensembl_transcript_id, .keep_all = TRUE) %>% dplyr::select(-start, -end)

# Find stop codon sequences
transcript_coding <- transcript_coding %>% mutate(coding_length = str_length(coding), stop_codons = substr(coding, coding_length - 2, coding_length)) %>% filter(stop_codons %in% c("TAA", "TGA", "TAG")) %>% dplyr::select(-coding)

# Process 3' UTR sequences
columns_3utr <- c("3_utr_start", "3_utr_end", "3utr")
rename_columns_3utr <- c("start_3utr", "end_3utr", "sequence_3utr")

# Function to translate DNA sequences
translate_sequence <- function(seq) {
  if (is.na(seq) || seq == "") {return("")}
  else {return(as.character(translate(DNAString(seq), no.init.codon = TRUE, if.fuzzy.codon = "error")))}
}

transcript_3utr <- transcript_3utr %>% rename_at(vars(columns_3utr), ~ rename_columns_3utr) %>% filter(sequence_3utr != "Sequence unavailable")
transcript_3utr <- transcript_3utr %>% mutate(length_3utr = str_length(sequence_3utr), gc_3utr = (str_count(sequence_3utr, "G") + str_count(sequence_3utr, "C")) / length_3utr * 100,
                                              translate_tail = sapply(sequence_3utr, translate_sequence),  translate_tail = ifelse(grepl("\\*", translate_tail), sapply(strsplit(translate_tail, "\\*"), function(x) paste0(x[1], "*")), NA),
                                              ext_length = ifelse(is.na(translate_tail), NA, nchar(translate_tail)), next_stop_codon = substr(sequence_3utr, ext_length*3-2, ext_length*3),
                                              shortened_3utr_ratio = (length_3utr - ext_length*3)/length_3utr * 100) %>% dplyr::select(-sequence_3utr)
  

# Combine coding and UTR information
sequence_info <- transcript_coding %>% inner_join(transcript_3utr, by = "ensembl_transcript_id") %>% dplyr::select(-strand)
MANE_transcripts <- inner_join(MANE_transcripts, sequence_info, by = "ensembl_transcript_id") %>% arrange(chromosome, POS_1)

# Function to get complementary sequence for negative strand
complement <- function(sequence) {chartr("ATGC", "TACG", sequence)}

MANE_transcripts <- MANE_transcripts %>% mutate(stop_codon_sequences = if_else(strand == "-", complement(stop_codons), stop_codons),
    REF_1 = substr(stop_codon_sequences, 1, 1), REF_2 = substr(stop_codon_sequences, 2, 2), REF_3 = substr(stop_codon_sequences, 3, 3))

# Extract variant position and convert to a VCF file
vcf_data_list <- list(
  MANE_transcripts %>% dplyr::select(chromosome, POS_1, REF_1),
  MANE_transcripts %>% dplyr::select(chromosome, POS_2, REF_2),
  MANE_transcripts %>% dplyr::select(chromosome, POS_3, REF_3)
)

# Function to generate alternate nucleotides
mutagenesis <- function(ref) {nucleotides <- c("A", "T", "G", "C")
  alt <- nucleotides[nucleotides != ref]
  return(alt)
}

vcf_data <- bind_rows(lapply(vcf_data_list, function(df) {
  names(df) <- c("chromosome", "POS", "REF")
  df %>% rowwise() %>% mutate(ALT = list(mutagenesis(REF))) %>% unnest(cols = c(ALT))
})) %>% arrange(chromosome, POS) %>% unique()

# Write the final summary table and vcf_data to files
write.table(MANE_transcripts, "mane_transcripts_info.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(vcf_data, vcf_body, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
