#!/bin/bash

sample_name=$1 #e.g.AZ107
genome_name=$2	#reference genome to be used

# Split bin into contigs
		
while read first; do
	if [ ${first:0:1} == ">" ] ; then
    		filename=$(echo "$first" | cut -d ":" -f1 | tr -d ">")
    		touch ./"${genome_name}"_"$filename".fa
   		echo "$first" >> ./"${genome_name}"_"${filename}".fa
  	else
   		echo "$first" >> ./"${genome_name}"_"${filename}".fa
  	fi
done < "${genome_name}.fa"
		
# Align reads to contigs
bowtie2-build "${genome_name}.fa" "${genome_name}"
bowtie2 \
	-t \
	-x "${genome_name}" \
	-1 "${sample_name}_filtered_pair1.fq" \
	-2 "${sample_name}_filtered_pair2.fq" \
	-p 30 \
	-S "${genome_name}.sam"

# Output only aligned reads in bam format
samtools view -bSF4 "${genome_name}.sam" > "${genome_name}.bam" -@ 30

# Sort bam file
samtools sort "${genome_name}.bam" -o "${genome_name}.sorted.bam"

# Index BAM file
samtools index "${genome_name}.sorted.bam"

# Count number of reads aligned to each contig
samtools idxstats "${genome_name}.sorted.bam" > "${genome_name}_num_reads_per_contig.txt"

# Based on third column, remove rows that have # of reads <1000 
awk -F"\t" '$3>999' "${genome_name}_num_reads_per_contig.txt" | cut -f1 > "${genome_name}_contigs.txt"

while read contigs;
do
	contig=$(echo "$contigs")

	# create new header
	samtools view -H "${genome_name}.sorted.bam" | grep -P -v "^@SQ" > "${genome_name}"_"${contig}.header"
	samtools view -H "${genome_name}.sorted.bam" | grep -P "^@SQ\tSN:${contigs}\t" >> "${genome_name}"_"${contig}.header"
	samtools view -bh "${genome_name}.sorted.bam" "${contig}" > "${genome_name}"_"${contig}.bam"
	(cat "${genome_name}"_"${contig}.header" <(samtools view "${genome_name}"_"${contig}.bam") | samtools view -bo "${genome_name}"_"${contig}_rehead.bam" -)

	# run DamageProfiler per contig
	java -jar DamageProfiler-0.4.7.jar -i "${genome_name}"_"${contig}_rehead.bam" -o out_"${genome_name}" -r "${genome_name}"_"${contig}.fa" -title "${genome_name}"_"${contig}" -yaxis_damageplot 0.10 2> "${genome_name}"_"${contig}_error"
		
done < "${genome_name}_contigs.txt"