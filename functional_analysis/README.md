## Overview
Pipeline to predict genes and CAZymes and process gene catalog.

## Files
### prokka.sh
Script to predict genes from contigs.

### cd-hit.sh
Build gene catalog from predicted genes.

### align_and_build_abundance_matrix.sh
Align raw reads to gene catalog and build a matrix of gene relative abundances.

### build_prokka_binary_and_count_matrix.sh
Create a binary matrix and a count matrix from PROKKA output files. We used this to build matrices to map functions to taxonomies.

### cazy_hmmscan.sh
Script to predict CAZymes from PROKKA protein output files by running hmmsearch against dbCAN HMMs v8.

### build_cazy_abundance_matrix.sh
Build CAZyme relative abundance matrix by dividing the number of times each CAZy family is predicted in each sample by the total number of CAZymes predicted in the sample.

### jaccard_distance.py
Calculate pairwise Jaccard distances.
