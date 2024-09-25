### `extract_ensembl_data.R`

This R script is designed to extract specific data from the Ensembl database, focusing on stop-lost variants and associated transcript information.

#### Required R Packages:
To run this script, you need to have the following R packages installed:
- `tidyverse`, `biomaRt`, `stringr`, `Biostrings`

#### Output Files:
The script generates the following output files:
1. **`mane_transcripts_info.txt`**: Contains information on stop codons for all MANE select transcripts retrieved from the Ensembl database.
2. **`stoplost_SNV.txt`**: All possible single nucleotide substitutions at stop codons, including chromosome, position, reference, and alternative sequences.
