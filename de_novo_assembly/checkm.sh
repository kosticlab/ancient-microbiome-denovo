#!/bin/bash

input_folder=$1 #directory containing all bins (from MetaBAT2 output)
sample_name=$2 #e.g.AZ107

#run CheckM
checkm lineage_wf \
	-f checkm_output.txt \
	-t 30 \
	-x fa \
	"${input_folder}" checkm_output

checkm qa \
	-o 2 \
	-t 30 \
	-f checkm_stats \
	--tab_table checkm_output/lineage.ms \
	checkm_output

#calculate coverage per contig
checkm coverage \
	"${sample_name}_final.contigs.fa.metabat-bins0" \
	-x fa \
	-t 35 \
	"${sample_name}_coverage.tsv" \
	"${sample_name}.sorted.bam"


##calculate coverage per genome by averaging the coverage levels across all contigs per genome
#make a list of bins per sample
ls "${sample_name}_final.contigs.fa.metabat-bins0" > "${sample_name}_bin_names"

# per sample, loop through all the bins
while read line;
do
	bin_name=$(echo "$line")
        echo $bin_name

	# create a file with bin names
	echo "${sample_name}"_"${bin_name}" >> "${sample_name}_bin_names.txt"

	# extract contigs from each bin that are not unbinned and sort by coverage
        awk '(NR>1) && ($2 == "'${bin_name/.fa/}'") ' "${sample_name}_coverage.tsv" | sort -k5 > "${sample_name}"_"${bin_name/.fa/}_coverage.tsv"

	# find average coverage per bin
	awk -F'\t' '{ sum += $5 } END { print sum / NR }' "${sample_name}"_"${bin_name/.fa/}_coverage.tsv" >> "${sample_name}_avg_cov_per_bin.txt"

done < "${sample_name}_bin_names" #list of all the bins per sample	


# paste the 2 lists together
paste "${sample_name}_bin_names.txt" "${sample_name}_avg_cov_per_bin.txt" > "${sample_name}_coverage_per_bin.txt"