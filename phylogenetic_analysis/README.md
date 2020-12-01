## Overview
Pipeline for building phylogenetic trees.

## Files
### phylogenetic_tree.sh
This script was used to build phylogenetic trees. The steps are as follows:
1. Run GTDB-Tk "classify" workflow to identify 120 bacterial marker genes and build a multiple sequence alignment (MSA) based on these marker genes.
2. The resulting FASTA files containing MSA of the submitted genomes were used for maximum likelihood phylogenetic tree inference using IQ-TREE.
