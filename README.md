# TAILVAR (Terminal codon Analysis and Improved prediction of Lengthened VARiants)
😊 Welcome to the **TAILVAR** repository! This repository stores the codes for developing the TAILVAR score designed to assess the functional impact of **stop-loss variants** occurring at stop codons (TAA, TGA, TAG) 🚀


![TAILVAR overview](images/TAILVAR_overview.jpg)

# Overview
**TAILVAR** is built using a Random Forest model that predicts the pathogenicity of **stop-loss variants**. By integrating a combination of in-silico prediction scores, transcript, and protein features of C-terminal extensions, **TAILVAR** provides a score ranging from 0 to 1, indicating the probability of pathogenic potential. TAILVAR score cutoffs ≥ 0.70 and ≤ 0.30 can be used to distinguish potential pathogenic/likely pathogenic and benign/likely benign variants.

![TAILVAR overview](images/TAILVAR_performance.jpg)

## Key components (integration of 36 features)

- **Variant effect prediction tool**:
  - **[CADD](http://cadd.gs.washington.edu/)** (Combined Annotation Dependent Depletion)

- **Conservation scores**:
  - **[GERP](http://mendel.stanford.edu/SidowLab/downloads/gerp/)** (Genomic Evolutionary Rate Profiling)
  - **[phyloP100way](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/)** (Phylogenetic P-value across 100 vertebrates)
  - **[phastCons100way](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/)** (Phylogenetic Conserved Elements across 100 vertebrates)

- **Transcript features**:
  - **3'UTR_GC**: GC content of the 3' UTR
  - **3'UTR_length**: Length of the 3' UTR
  - **[pLI]([http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way](https://gnomad.broadinstitute.org/data#v4-constraint)/)**: probablity score of loss-of-function (LOF) intolerance from gnomAD
  - **[LOEUF](https://gnomad.broadinstitute.org/data#v4-constraint)**: upper boundary fraction of observed/expected LOF variants from gnomAD
  - **[mRNA stability](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02811-x)**: Z-scores of Saluki dataset

- **Protein features**:
  - **Protein_lengths**: Total counts of amino acids in the original protein
  - **C-terminal_lengths**: Total counts of amino acids in the C-terminal extension
  - **Amino acids counts**: Counts of each of the 20 amino acids in the C-terminal extension
  - **Hydrophobicity_KD**: Mean hydrophobicity of the C-terminal extension in Kyte–Doolittle (KD) scales
  - **Hydrophobicity_MJ**: Mean hydrophobicity of the C-terminal extension in Miyazawa–Jernigan (MJ) scales
  - **[TANGO](https://tango.crg.es/)**: aggregation properties of the C-terminal peptide predicted by TANGO
  - **[CANYA](https://github.com/lehner-lab/canya)**: aggregation properties of the C-terminal peptide predicted by CANYA   

# Annotation
To annotate **TAILVAR** scores in **[VEP](https://github.com/Ensembl/ensembl-vep)**, use the following command:
```bash
./vep [...] --custom file=TAILVAR_score.vcf.gz,short_name=TAILVAR,format=vcf,type=exact,fields=score
```

# Download
The annotated dataset including **TAILVAR** scores for all possible single nucleotide substitutions are available for download [here](https://zenodo.org/records/14512203).

# Citation
If you use **TAILVAR** in your research, please cite:

JG Yoon. *A predictive framework for stop-loss variants with C-terminal extensions.*
