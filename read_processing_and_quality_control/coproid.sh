#!/bin/bash

nextflow run nf-core/coproid \
	--reads '*paired_{1,2}.fastq.gz' \
	--name1 "Homo_sapiens" \
	--name2 "Canis_familiaris" \
	--genome1 'Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa' \
	--genome2 'Canis_familiaris/NCBI/build3.1/Sequence/WholeGenomeFasta/genome.fa' \
	--krakendb '/home/ubuntu/' \
	--outdir coproid_output \
	--max_cpus 15 \
	-profile conda \
	-c ~/.nextflow/config
