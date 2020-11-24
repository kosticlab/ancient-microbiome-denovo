#!/bin/bash

# Run the first part of the script for each sample

sample_name=$1 #e.g.AZ107


### Process hmmsearch results
# take the first column (target name) and fifth column (e-value for full sequence)
tail -n +4 "${sample_name}_cazymes_table.txt" | head -n -10 | sed 's/ \{2,\}/\t/g' | cut -f1,3,5 > "${sample_name}_cazy_groups.txt"

# remove rows with e-value >=1e-5
awk '$3 <1e-5' "${sample_name}_cazy_groups.txt" > "${sample_name}_cazy_groups_cutoff.txt"

# add header
cat cazyme_header.txt "${sample_name}_cazy_groups_cutoff.txt" > "${sample_name}_cazy_groups_cutoff_header.tsv"


### Create an abundance matrix

# get just the protein names per sample
cut -f1,2 "${sample_name}_cazy_groups_cutoff_header.tsv" > foo1

# count unique cazy group names
cut -f2 foo1 | tail -n +2 | sort | uniq -c > test_compiled3.txt

# Replace spaces one or more with a tab
sed -i 's/ \{1,\}/\t/g' test_compiled3.txt
cut -f2 test_compiled3.txt > counts
cut -f3 test_compiled3.txt > prot
paste prot counts > "${sample_name}_prot_counts"



### Run this second part for all of the samples together

# merge tables using merge_metaphlan_tables.py
python merge_metaphlan_tables.py *prot_counts > cazy_groups_merged.txt
