#!/bin/bash

# run this script for each sample

sample_name=$1 #e.g.AZ107

# run PROKKA
/mnt/pkgs/prokka/bin/prokka \
	--notbl2asn \
	--addgenes \
	--metagenome \
	--cpus 0 \
	--mincontiglen 1 \
	--outdir "${sample_name}_output"/prokka-output \
	"${sample_name}_output"/megahit-output/final.contigs.fa \
	2>"${sample_name}_output"/"${ID1}_error.log"

