#!/bin/bash

#run this script per sample

sample_name=${1} #e.g.AZ107

#run AdapterRemoval
AdapterRemoval --file1 "${sample_name}_1.fq" --file2 "${sample_name}_2.fq" --basename "${sample_name}" --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --trimns --trimqualities --threads 15

#run KneadData
kneaddata --input "${sample_name}_1.fq" --input /home/ubuntu/"${sample_name}_2.fq" --reference-db /mnt/pkgs/knead_database/ --threads 30 --output "${sample_name}" --bypass-trim

#run Cutadapt
cutadapt --cores=35 --minimum-length 30 -o "${sample_name}_filtered_pair1.fq" -p "${sample_name}_filtered_pair2.fq" "${sample_name}_kneaddata_paired_1.fastq" "${sample_name}_kneaddata_paired_2.fastq" &>"${sample_name}_out_messages.txt"

