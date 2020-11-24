#!/bin/bash

input_folder=$1 #directory containing all MAGs to be placed on phylogenetic tree

#run CheckM
gtdbtk classify_wf \
	--cpus 70 \
	--extension fa \
	--genome_dir "${input_folder}" \
	--out_dir gtdbtk_output

iqtree -s gtdbtk_output/gtdbtk.bac120.user_msa.fasta -nt AUTO -m LG
