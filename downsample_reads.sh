#!/bin/bash

sample_name=$1 #e.g.AZ107
num_reads=$2 #desired final number of reads

seqtk sample -s5 "${sample_name}_filtered_pair1.fq" "${num_reads}" > "${sample_name}_downsampled_1.fq"
seqtk sample -s5 "${sample_name}_filtered_pair2.fq" "${num_reads}" > "${sample_name}_downsampled_2.fq"