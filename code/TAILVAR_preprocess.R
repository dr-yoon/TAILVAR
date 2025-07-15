# Load necessary libraries
library(tidyverse)
library(Peptides)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]

# Load transcript info data
transcript_info <- read_tsv("mane_transcripts_info.txt")
peptide_info <- read_tsv("mane_transcripts_peptide.txt")
transcript_info <- transcript_info %>% left_join(peptide_info, by = c("ensembl_transcript_id")) %>% distinct()
transcript_info <- transcript_info %>% dplyr::select(ensembl_transcript_id, length_3utr, gc_3utr, stop_codons, next_stop_codon, coding_length, peptide, translate_tail)
transcript_info_colname <- c("ensembl_transcript_id", "length_3utr", "gc_3utr", "coding_length", "peptide", "translate_tail")
transcript_info_rename <- c("Feature", "UTR3_length", "UTR3_GC","protein_length", "protein_sequence", "extension_peptide")

transcript_info <- transcript_info %>% rename_at(vars(transcript_info_colname), ~ transcript_info_rename)

# Load dataset
data_input <- read_tsv(paste0(file_name,"_filtered.txt"))
data_input <- left_join(data_input, transcript_info, by = c("Feature")) %>% distinct()
data_input <- data_input %>% mutate(extension_peptide = paste0(Amino_acids, extension_peptide), extension_peptide = str_remove_all(extension_peptide, fixed("*")))
data_input <- data_input %>% mutate(extended_length = as.numeric(extended_length)) %>% filter(is.na(extended_length) == FALSE)

gnomad <- read_tsv("gnomad.v4.1.constraint_metrics.tsv") %>% dplyr::select(transcript, lof.pLI, lof.oe_ci.upper)
saluki <- read_tsv("Saluki_data.txt") %>% dplyr::select(gene_id, Saluki)
IDP_data <- read_tsv("DisProt_Human_IDP_info.txt")

data_input <- data_input %>% left_join(gnomad, by = c("Feature" = "transcript")) %>% rename_with(~ c("Clinvar", "pLI","LOEUF"), c("ClinVar_CLNSIG","lof.pLI","lof.oe_ci.upper"))
data_input <- data_input %>% left_join(saluki, by = c("Gene" = "gene_id")) %>% rename_with(~ c("mRNA_stability"), c("Saluki"))
data_input <- data_input %>% mutate(IDP = if_else(Uniprot %in% IDP_data$acc, 1, 0))

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
deeploc_data <- data_input %>% dplyr::select("SYMBOL", "protein_sequence", "Var_ID", "extension_peptide") %>% distinct()
deeploc_data <- deeploc_data %>% mutate(WT = gsub("\\*$", "", protein_sequence), MT = paste0(WT,extension_peptide))
fasta_deeploc <- apply(deeploc_data, 1, function(row) {
  wt_header <- paste0(">", row["SYMBOL"], "_WT")
  mt_header <- paste0(">", row["Var_ID"])
  wt_seq <- row["WT"]
  mt_seq <- row["MT"]
  paste(c(wt_header, wt_seq, mt_header, mt_seq), collapse = "\n")
})
writeLines(fasta_deeploc, paste0(file_name,".deeploc.fasta"))

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
