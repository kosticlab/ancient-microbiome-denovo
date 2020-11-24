#!/bin/bash

#run this script per sample

sample_name=$1 #e.g.AZ107

/mnt/pkgs/metaphlan2/metaphlan2.py \
	"${sample_name}_filtered_pair1.fq","${sample_name}_filtered_pair2.fq" \
	-o "${sample_name}_filtered_metaphlan" \
	--input_type fastq \
	--bowtie2out "${sample_name}_filtered_bowtieoutput" \
	--nproc 35
