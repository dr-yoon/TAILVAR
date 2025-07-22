# Load necessary libraries
library(tidyverse)
library(biomaRt)
library(stringr)
library(Biostrings)
library(Peptides)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]

# Load dataset
MANE_transcripts <- read_tsv("MANE_stop_codon_info.tsv") %>% mutate(ensembl_transcript_id = sub("\\..*$", "", transcript_id)) %>% dplyr::select(-transcript_id)
data_input <- read_tsv(paste0(file_name,"_filtered.txt"))

gnomad <- read_tsv("gnomad.v4.1.constraint_metrics.tsv") %>% dplyr::select(transcript, lof.pLI, lof.oe_ci.upper)
saluki <- read_tsv("Saluki_data.txt") %>% dplyr::select(gene_id, Saluki)
IDP_data <- read_tsv("DisProt_Human_IDP_info.txt")

data_input <- data_input %>% left_join(gnomad, by = c("Feature" = "transcript")) %>% rename_with(~ c("Clinvar", "pLI","LOEUF"), c("ClinVar_CLNSIG","lof.pLI","lof.oe_ci.upper"))
data_input <- data_input %>% left_join(saluki, by = c("Gene" = "gene_id")) %>% rename_with(~ c("mRNA_stability"), c("Saluki"))
data_input <- data_input %>% mutate(IDP = if_else(Uniprot %in% IDP_data$acc, 1, 0))

# Obtain coding and 3'UTR data from the Ensembl database
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#transcript_coding <- getBM(attributes = c("ensembl_transcript_id", "coding"), filters = "ensembl_transcript_id", values = MANE_transcripts$ensembl_transcript_id, mart = ensembl)
#transcript_3utr <- getBM(attributes = c("ensembl_transcript_id", "3_utr_start", "3_utr_end", "strand", "3utr"), filters = "ensembl_transcript_id", values = MANE_transcripts$ensembl_transcript_id, mart = ensembl)

# Read pre-downloaded data: since querying the Ensembl database can be time-consuming, we use pre-downloaded data files
transcript_coding <- read_tsv("MANE_transcripts_coding.txt") %>% dplyr::select("ensembl_transcript_id", "coding") %>%
                  mutate(protein_length = nchar(coding)/3 - 1, stop_codons = substr(coding, nchar(coding)-2, nchar(coding))) %>% filter(stop_codons %in% c("TGA", "TAA", "TAG"))
transcript_3utr <- read_tsv("MANE_transcripts_3utr.txt") %>% dplyr::select("ensembl_transcript_id", "3utr") %>% dplyr::rename(sequence_3utr = `3utr`) %>% filter(sequence_3utr != "Sequence unavailable")
transcript_3utr <- transcript_3utr %>% mutate(UTR3_length = str_length(sequence_3utr), UTR3_GC = (str_count(sequence_3utr, "G") + str_count(sequence_3utr, "C")) / UTR3_length * 100)
       
data_input <- inner_join(data_input, transcript_coding, by = c("Feature" = "ensembl_transcript_id"))
data_input <- inner_join(data_input, transcript_3utr, by = c("Feature" = "ensembl_transcript_id"))

# Function to translate C-terminal pepetide sequences
translate_sequence <- function(seq) {
  if (is.na(seq) || seq == "") return("")
  pep <- as.character(translate(DNAString(seq), no.init.codon = TRUE, if.fuzzy.codon = "error"))
  strsplit(pep, "\\*", fixed = FALSE)[[1]][1]
}

data_input <- data_input %>% mutate(extended_ORF = paste0(Codons_ALT,sequence_3utr)) %>%
  mutate (extension_peptide = vapply(extended_ORF, translate_sequence, character(1)), Amino_acids = substr(extension_peptide, 1, 1),
          extended_length = ifelse(is.na(extension_peptide), NA, nchar(extension_peptide)),
          next_stop_codon = substr(extended_ORF, extended_length*3+1, extended_length*3+3),
          next_stop_codon = if_else(!next_stop_codon %in% c("TAG", "TGA", "TAA"), NA, next_stop_codon),
          extended_length = if_else(is.na(next_stop_codon), NA, extended_length)) %>% dplyr::select(-sequence_3utr, -coding, -extended_ORF)

data_input <- data_input %>% mutate(extended_length = as.numeric(extended_length)) %>% filter(!is.na(extended_length) & extended_length != 0)


# Add columns for each amino acid count
amino_acid_columns <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (aa in amino_acid_columns) {data_input <- data_input %>% mutate(!!aa := str_count(extension_peptide, aa))}

hydro_scale <- function(df) {
  df %>% rowwise() %>% 
    mutate(
      hydro_KD = if (is.na(extension_peptide) || extension_peptide == "") NA_real_
      else hydrophobicity(extension_peptide, scale = "KyteDoolittle"),
      hydro_MJ = if (is.na(extension_peptide) || extension_peptide == "") NA_real_
      else hydrophobicity(extension_peptide, scale = "Miyazawa")
    ) %>% ungroup()
}

data_input <- hydro_scale(data_input) %>% mutate(Var_ID = paste0(SYMBOL,"_",HGVSp))

# Save peptide sequence input
seq_data <- data_input %>% dplyr::select("Var_ID","extension_peptide") %>% distinct()
write.table(seq_data, paste0(file_name,"_seq_input.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# Define features and target
Pop_AF <-c("RegeneronME_ALL_AF","AllofUs_POPMAX_AF","ALFA_Total_AF","gnomAD4.1_AF_joint")
comp_scores <- c("CADD", "DANN", "Eigen_PC", "BayesDel_noAF","BayesDel_addAF","FATHMM_MKL","fitCons_int", "MutationTaster","VEST4","GPN_MSA", "GERP", "phyloP100", "phastCons100")

Model_variables  <- c("Upload_variation", "Feature", "SYMBOL", "HGVSc", "HGVSp", "pLI", "LOEUF", all_of(Pop_AF), all_of(comp_scores),
                      "Clinvar", "UTR3_length", "UTR3_GC", "protein_length", "extended_length", "extension_peptide", "stop_codons", "next_stop_codon",
                      all_of(amino_acid_columns), "Amino_acids", "hydro_KD", "hydro_MJ","Uniprot", "mRNA_stability","IDP","Var_ID")

preprocessed_data <- data_input %>% dplyr::select(all_of(Model_variables))

# Write the filtered dataset to a file
write.table(preprocessed_data, paste0(file_name,"_preprocessed.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
