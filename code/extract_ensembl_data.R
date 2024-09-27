# Load necessary libraries
library(tidyverse)
library(biomaRt)
library(stringr)
library(Biostrings)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define attributes for querying the MANE project
attributes_list <- list(
  c("ensembl_transcript_id", "transcript_mane_select", "external_gene_name", "chromosome_name", "strand", "percentage_gene_gc_content"),
  c("ensembl_transcript_id", "chromosome_name", "genomic_coding_start", "genomic_coding_end", "strand", "rank"),
  c("ensembl_transcript_id", "coding"),
  c("ensembl_transcript_id", "5utr"),
  c("ensembl_transcript_id", "3utr")
)

# Query the Ensembl database download
transcript_info_1 <- getBM(attributes = attributes_list[[1]], filters = "mane_select", values = TRUE, mart = ensembl)
transcript_info_2 <- getBM(attributes = attributes_list[[2]], filters = "mane_select", values = TRUE, mart = ensembl)
transcript_info_3 <- getBM(attributes = attributes_list[[3]], filters = "mane_select", values = TRUE, mart = ensembl)
transcript_info_4 <- getBM(attributes = attributes_list[[4]], filters = "mane_select", values = TRUE, mart = ensembl)
transcript_info_5 <- getBM(attributes = attributes_list[[5]], filters = "mane_select", values = TRUE, mart = ensembl)

# Read pre-downloaded data: since querying the Ensembl database can be time-consuming, 
# we use pre-downloaded data files instead. Uncomment the lines below to read the files.
# transcript_info_1 <- read_tsv("mane_transcripts_gene.txt")
# transcript_info_2 <- read_tsv("mane_transcripts_position.txt")
# transcript_info_3 <- read_tsv("mane_transcripts_coding.txt")
# transcript_info_4 <- read_tsv("mane_transcripts_5utr.txt")
# transcript_info_5 <- read_tsv("mane_transcripts_3utr.txt")

# Filter for valid chromosomes and remove NA values
valid_chromosomes <- c(as.character(1:22), "X", "Y")
transcript_info_1 <- transcript_info_1 %>% filter(chromosome_name %in% valid_chromosomes)
transcript_info_2 <- transcript_info_2 %>% filter(chromosome_name %in% valid_chromosomes, !is.na(genomic_coding_end)) %>%
  mutate(coding_length = genomic_coding_end - genomic_coding_start + 1)

# Identify splitted transcripts
splitted_transcript <- transcript_info_2 %>%
  group_by(ensembl_transcript_id) %>%
  top_n(1, wt = rank) %>%
  filter(coding_length <= 1)

# Identify and summarize non-splitted transcripts
no_splitted_transcript <- transcript_info_2 %>%
  filter(!ensembl_transcript_id %in% splitted_transcript$ensembl_transcript_id) %>%
  group_by(ensembl_transcript_id) %>% summarize(
    chromosome = dplyr::first(chromosome_name), strand = dplyr::first(strand),
    coding_start = min(genomic_coding_start, na.rm = TRUE), coding_end = max(genomic_coding_end, na.rm = TRUE)
  ) %>% mutate(
    POS_1 = if_else(strand == 1, coding_end - 2, coding_start),
    POS_2 = if_else(strand == 1, coding_end - 1, coding_start + 1),
    POS_3 = if_else(strand == 1, coding_end, coding_start + 2)
  ) %>% dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Identify and summarize splitted transcripts
splitted_transcript2 <- transcript_info_2 %>%
  filter(ensembl_transcript_id %in% splitted_transcript$ensembl_transcript_id) %>%
  group_by(ensembl_transcript_id) %>% top_n(2, wt = rank) %>%
  arrange(ensembl_transcript_id, desc(rank)) %>%
  summarize(
    chromosome = dplyr::first(chromosome_name), strand = dplyr::first(strand),
    last_exon_overlap = genomic_coding_end[rank == max(rank)] - genomic_coding_start[rank == max(rank)] + 1,
    exon_junction_3end = if_else(strand == 1, genomic_coding_end[rank == min(rank)], genomic_coding_end[rank == max(rank)]),
    exon_junction_5end = if_else(strand == 1, genomic_coding_start[rank == max(rank)], genomic_coding_start[rank == min(rank)])
  ) %>% mutate(
    POS_1 = if_else(strand == 1, if_else(last_exon_overlap == 1, exon_junction_3end - 1, exon_junction_3end), if_else(last_exon_overlap == 1, exon_junction_3end, exon_junction_3end - 1)),
    POS_2 = if_else(strand == 1, if_else(last_exon_overlap == 1, exon_junction_3end, exon_junction_5end), if_else(last_exon_overlap == 1, exon_junction_5end, exon_junction_3end)),
    POS_3 = if_else(strand == 1, if_else(last_exon_overlap == 1, exon_junction_5end, exon_junction_5end + 1), if_else(last_exon_overlap == 1, exon_junction_5end + 1, exon_junction_5end))
  ) %>% dplyr::select(ensembl_transcript_id, POS_1, POS_2, POS_3)

# Combine non-splitted and splitted transcript positions
stop_codon_pos <- bind_rows(no_splitted_transcript, splitted_transcript2)

# Merge transcript information with stop codon positions
transcript_info <- left_join(transcript_info_1, stop_codon_pos, by = "ensembl_transcript_id")

# Function to translate DNA sequences
translate_sequence <- function(seq) {
  if (is.na(seq) || seq == "") {return("")}
  else {return(as.character(translate(DNAString(seq), no.init.codon = TRUE, if.fuzzy.codon = "error")))}
}

# Calculate coding and UTR lengths and GC contents
transcript_info_3 <- transcript_info_3 %>% mutate(
    coding_length = str_length(coding), 
    gc_coding = (str_count(coding, "G") + str_count(coding, "C")) / coding_length * 100,
    stop_codons = substr(coding, coding_length - 2, coding_length),
    peptide = sapply(coding, translate_sequence),
  )

# Process 5' UTR sequences
names(transcript_info_4)[2] <- "sequence_5utr"
transcript_info_4 <- transcript_info_4 %>% filter(sequence_5utr != "Sequence unavailable") %>%
  mutate(
    length_5utr = str_length(sequence_5utr), 
    gc_5utr = (str_count(sequence_5utr, "G") + str_count(sequence_5utr, "C")) / length_5utr * 100
  )

# Process 3' UTR sequences
names(transcript_info_5)[2] <- "sequence_3utr"
transcript_info_5 <- transcript_info_5 %>% filter(sequence_3utr != "Sequence unavailable") %>%
  mutate(
    length_3utr = str_length(sequence_3utr), 
    gc_3utr = (str_count(sequence_3utr, "G") + str_count(sequence_3utr, "C")) / length_3utr * 100
  )

transcript_info_5 <- transcript_info_5 %>%
  mutate(translate_tail = sapply(sequence_3utr, translate_sequence),
    translate_tail = ifelse(grepl("\\*", translate_tail), sapply(strsplit(translate_tail, "\\*"), function(x) paste0(x[1], "*")), NA),
    ext_length = ifelse(is.na(translate_tail), NA, nchar(translate_tail)), next_stop_codon = substr(sequence_3utr, ext_length*3-2, ext_length*3))

# Combine coding and UTR information
sequence_info <- transcript_info_3 %>% left_join(transcript_info_4, by = "ensembl_transcript_id") %>% left_join(transcript_info_5, by = "ensembl_transcript_id")
transcript_info <- left_join(transcript_info, sequence_info, by = "ensembl_transcript_id") %>% arrange(chromosome_name, POS_1)

# Function to get reverse complementary sequence
reverse_complement <- function(sequence) {chartr("ATGC", "TACG", sapply(strsplit(sequence, NULL), function(x) paste(rev(x), collapse = "")))}

transcript_info <- transcript_info %>% mutate(
    stop_codon_sequences = if_else(strand == -1, reverse_complement(stop_codons), stop_codons),
    REF_1 = substr(stop_codon_sequences, 1, 1),
    REF_2 = substr(stop_codon_sequences, 2, 2),
    REF_3 = substr(stop_codon_sequences, 3, 3)
  )

# Extract variant position and convert to a VCF file
vcf_data_list <- list(
  transcript_info %>% dplyr::select(chromosome_name, POS_1, REF_1),
  transcript_info %>% dplyr::select(chromosome_name, POS_2, REF_2),
  transcript_info %>% dplyr::select(chromosome_name, POS_3, REF_3)
)

# Function to generate alternate nucleotides
mutagenesis <- function(ref) {nucleotides <- c("A", "T", "C", "G")
  alt <- nucleotides[nucleotides != ref]
  return(alt)
}

vcf_data <- bind_rows(lapply(vcf_data_list, function(df) {
  names(df) <- c("chromosome", "POS", "REF")
  df %>% rowwise() %>% mutate(ALT = list(mutagenesis(REF))) %>% unnest(cols = c(ALT))
})) %>% arrange(chromosome, POS) %>% mutate(chromosome = paste0("chr", chromosome))

# Write the final summary table and vcf_data to files
write.table(transcript_info, "mane_transcripts_info.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(vcf_data, "stoplost_SNV.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
