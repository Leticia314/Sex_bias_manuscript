# Sex bias project

This repository contains the code used in the manuscript "Sex-biased gene expression across mammalian organ development and evolution".

- The raw data generated in this study are deposited in ArrayExpress with the accession code E-MTAB-12180.

- Processed data missing in this repository (due to space constrains) can be downloaded from: https://heibox.uni-heidelberg.de/d/34ace5842ae34f798277/.

- Sex-biased genes can be explored interactively in our shiny app: https://apps.kaessmannlab.org/sexbiasapp.


# Data folders

- Raw_counts contains tables with raw counts for each species.

- Norm_counts contains tables with normalized counts for each species.

- RPKMs contains tables with rpkms for each species.

- Matadata contains tables with metadata information for each species.

- Calls_all_methods contains tables with how many and which tools called genes as "sex-biased" for each species.

- SB_genes_tables contains tables with all information regarding sex bias for each gene in each species.


# Code folders

- Simulations contains code regarding section "Simulation of sex-biased dataset" from the manuscript.

- Calling_SB_genes contains code for running all four tools for detecting sex-biased genes across development.

- DESeq2_adults contains code for detecting sex-biased genes in adult stages.

- Onset_SB_expression contains code for clustering sex-biased genes according to temporal dynamics.

- Single_cell_analyses contains example code for looking at the expression of sex-biased genes in sc/snRNAseq datasets.

- ChIP_seq_analyses contains example code regarding the section "Moluecular basis of sex-biased gene expression".
