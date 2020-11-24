#!/bin/bash

input_folder=$1 #directory containing all ancient MAGs (from MetaBAT2 output)
reference_genomes_folder=$2	#directory containing all reference genomes
bin_names=$3	#file listing the names of all the ancient MAGs
ancient_ref_bins=$4	#directory with all representative ancient MAGs and their top 100 closest reference genomes

## determine whether each SGB belong to a known microbial species
# Mash sketch for all of the reference genomes
mash sketch \
	-p 30 \
	-o mash_reference \
	"${reference_genomes_folder}"/*fa

# calculate Mash distance per MAG
while read line;
do
        name_of_bin=$(echo "$line")
        echo $name_of_bin
	
	mash dist -d 0.10 mash_reference.msh "${input_folder}"/"${name_of_bin}" | sort -gk3 | head -100 >> distance_top100_0.10.txt

done < "${bin_names}"


# calculate ANIs for each ancient MAG and its 100 closest reference genomes within 10% Mash distance
dRep cluster \
	-p 20 \
	-sa 0.95 \
	--S_algorithm fastANI \
	--SkipMash \
	-g "${ancient_ref_bins}"/*.fa drep_workdir 
