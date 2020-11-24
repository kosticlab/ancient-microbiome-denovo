#!/bin/bash

#run this script for each sample

# concatenate PROKKA output files from all samples
cat *ffn > all_prokka_combined.ffn


# run CD-HIT across all samples
cd-hit-est \
	-n 10 \
	-i all_prokka_combined.ffn \
	-o all_gene_catalos \
	-c .95 \
	-M 300000 \
	-s .9 \
	-T 0 \
	-aS .9

