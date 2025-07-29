# Load necessary libraries
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]

# Read the dataset
dataset <- read_tsv(paste0(file_name,"_anno.txt"), comment = "##")
transcript_info <- read_tsv("mane_transcripts_info.txt")

# Filter and process the dataset
names(dataset)[1] <- "Upload_variation"
dataset <- dataset %>% filter(Feature %in% transcript_info$ensembl_transcript_id, Codons != "-", str_detect(as.character(Consequence), "stop_lost")) %>%
  distinct() %>%  mutate(HGVSc = str_extract(HGVSc, "(?<=:).*"), HGVSp = str_extract(HGVSp, "(?<=:).*"), HGVSp = gsub("%3D", "=", HGVSp), Amino_acids = str_extract(Amino_acids, "(?<=/).*"), extended_length = str_extract(HGVSp, "(?<=fsTer|extTer)(\\d+|\\?)")) %>%
  separate(Codons, into = c("Codons_REF", "Codons_ALT"), sep = "\\/") %>% separate(Location, into = c("Chromosome", "Position"), sep = ":")

# Function to extract the max score
extract_score_max <- function(score_str) {
  if (is.na(score_str) | score_str == "-") { return(NA) }
  values <- str_split(score_str, ",", simplify = TRUE)
  values[values == "."] <- NA
  numeric_values <- as.numeric(values)
  if (all(is.na(numeric_values))) { return(NA) }
  max_values <- max(numeric_values, na.rm = TRUE)
  return(round(max_values, 3))
}

dataset <- dataset %>% mutate(MutationTaster_score = map_dbl(MutationTaster_score, extract_score_max), VEST4_score = map_dbl(VEST4_score, extract_score_max))

Pop_AF <-c("RegeneronME_ALL_AF","AllofUs_POPMAX_AF","ALFA_Total_AF","gnomAD4.1_nhomalt_joint", "gnomAD4.1_AF_grpmax_joint","gnomAD4.1_AF_joint")
comp_scores <- c("CADD_phred", "DANN_score", "Eigen-phred_coding", "BayesDel_noAF_score","BayesDel_addAF_score","fathmm-MKL_coding_score", "integrated_fitCons_score","MutationTaster_score","VEST4_score","GPN_MSA_score", "GERP++_RS", "phyloP100way_vertebrate","phastCons100way_vertebrate")

dataset <- dataset %>% 
  mutate(across(all_of(Pop_AF), ~as.numeric(.))) %>% mutate(across(all_of(Pop_AF), ~as.numeric(replace_na(., 0)))) %>%
  mutate(across(all_of(comp_scores), ~as.numeric(.)))
dataset$Uniprot <- sub("\\..*", "", dataset$SWISSPROT)

cols_remove <- c("Feature_type", "Existing_variation", "IMPACT", "DISTANCE", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "MANE_PLUS_CLINICAL", "TREMBL", "UNIPARC", "SWISSPROT","UNIPROT_ISOFORM", "INTRON", "HGVS_OFFSET", "ClinVar", "GPN_MSA", "gnomAD4.1")
dataset <- dataset %>% select(-any_of(cols_remove)) %>% arrange(Chromosome)

comp_scores_rename <- c("CADD", "DANN", "Eigen_PC", "BayesDel_noAF", "BayesDel_addAF", "FATHMM_MKL", "fitCons_int", "MutationTaster", "VEST4", "GPN_MSA", "GERP", "phyloP100", "phastCons100")
dataset <- dataset %>% rename_at(vars(comp_scores), ~ comp_scores_rename)

# Write the filtered dataset to a file
write.table(dataset, paste0(file_name,"_filtered.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
