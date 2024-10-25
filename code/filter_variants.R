# Load necessary libraries
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]

# Read the dataset
dataset <- read_tsv(paste0(file_name,"_anno.txt"), comment = "##")
transcript_info <- read_tsv("mane_transcripts_info.txt")

# Filter and process the dataset
names(dataset)[1] <- "Upload_variation"
dataset <- dataset %>% filter(Feature %in% transcript_info$ensembl_transcript_id, Codons != "-", str_detect(Consequence, "stop_lost")) %>%
  distinct() %>%  mutate(HGVSc = str_extract(HGVSc, "(?<=:).*"), HGVSp = str_extract(HGVSp, "(?<=:).*"), HGVSp = gsub("%3D", "=", HGVSp),
                         Amino_acids = str_extract(Amino_acids, "(?<=/).*"), extended_length = str_extract(HGVSp, "(?<=fsTer|extTer)(\\d+|\\?)")) %>%
  separate(Codons, into = c("Codons_REF", "Codons_ALT"), sep = "\\/") %>% separate(Location, into = c("Chromosome", "Position"), sep = ":")

Pop_AF <-c("gnomAD_exomes_AF", "gnomAD_exomes_POPMAX_AF", "gnomAD_genomes_AF", "gnomAD_genomes_POPMAX_AF","ALFA_European_AF","ALFA_African_AF","ALFA_Asian_AF","ALFA_Other_AF","ALFA_Total_AF")
comp_scores <- c("CADD_phred", "DANN_score", "fathmm-MKL_coding_score", "Eigen-phred_coding", "BayesDel_addAF_score", "BayesDel_noAF_score", "integrated_fitCons_score", "GERP++_RS", "phyloP100way_vertebrate", "phastCons100way_vertebrate")

dataset <- dataset %>% 
  mutate(across(all_of(Pop_AF), ~as.numeric(.))) %>% mutate(across(all_of(Pop_AF), ~as.numeric(replace_na(., 0)))) %>%
  mutate(across(all_of(comp_scores), ~as.numeric(.)))

dataset <- dataset %>% dplyr::select(-Gene, -Feature_type, -Existing_variation, -IMPACT, -DISTANCE, -FLAGS, -VARIANT_CLASS, -SYMBOL_SOURCE, -HGNC_ID, -BIOTYPE, -MANE_PLUS_CLINICAL, -SOURCE, -INTRON, -HGVS_OFFSET) %>% arrange(Chromosome)

comp_scores_rename <- c("CADD", "DANN", "FATHMM", "Eigen", "BayesDel_addAF", "BayesDel_noAF", "integrated_fitCons", "GERP", "phyloP100way", "phastCons100way")

dataset <- dataset %>% rename_at(vars(comp_scores), ~ comp_scores_rename)

# Write the filtered dataset to a file
write.table(dataset, paste0(file_name,"_filtered.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
