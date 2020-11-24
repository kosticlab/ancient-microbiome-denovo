#!/bin/bash

# run this script for each sample

sample_name=$1 #e.g.AZ107

hmmsearch \
	-o "${sample_name}_cazymes.txt" \
	--tblout "${sample_name}_cazymes_table.txt" \
	--cpu 20 \
	--incE 1e-3 \
	dbCAN-db/dbCAN-HMMdb-V8.txt \
	"${sample_name}_PROKKA_.faa"