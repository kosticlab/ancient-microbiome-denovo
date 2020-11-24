#!/bin/bash

#run this script for all samples

sample_name=$1 #e.g.AZ107

#run MEGAHIT
mkdir -p "${sample_name}_output"

/mnt/pkgs/megahit/megahit \
	-m 0.9 \
	-t 35 \
	-1 "${sample_name}_filtered_pair1.fq" \
	-2 "${sample_name}_filtered_pair2.fq" \
	-o "${sample_name}_output"/megahit-output

#build bowtie index
bowtie2-build "${sample_name}_output"/megahit-output/final.contigs.fa "${sample_name}"

#run bowtie2
bowtie2 \
	-t \
	-x "${sample_name}" \
	-1 "${sample_name}_filtered_pair1.fq" \
	-2 "${sample_name}_filtered_pair2.fq" \
	-p 30 \
	-S "${sample_name}.sam"

# Convert SAM to BAM file
samtools view -S -b "${sample_name}.sam" > "${sample_name}.bam" -@ 30

# Sort BAM file
samtools sort "${sample_name}.bam" -o "${sample_name}.sorted.bam"

# Index BAM file
samtools index "${sample_name}.sorted.bam"

# run MetaBAT2
runMetaBat.sh \
	-t 0 \
	"${sample_name}_output"/megahit-output/final.contigs.fa \
	"${sample_name}.sorted.bam"