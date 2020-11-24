#!/bin/bash

input_folder=$1 #directory containing all bins (from MetaBAT2 output)

# assign taxonomies using GTDB-Tk
gtdbtk classify_wf \
	--cpus 70 \
	--extension fa \
	--genome_dir "${input_folder}" \
	--out_dir gtdbtk_output_all
