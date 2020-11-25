#!/bin/bash

reference=${1} #reference genome
sample=${2}	#e.g.AZ107

# build bowtie index
bowtie2-build "${reference}.fa" "${reference}"

# run bowtie2
bowtie2 \
	-t \
	-x "${reference}" \
	-1 "${sample}_1.fq" \
	-2 "${sample}_2.fq" \
	-p 30 \
	-S "${sample}_${reference}.sam"

# run mapDamage
/home/ubuntu/mapDamage/bin/mapDamage \
	-y 0.04 \
	-t "${sample}_${reference}" \
	-i "${sample}_${reference}.sam" \
	-r "${reference}.fa" \
	2>"${sample}_${reference}_error"
