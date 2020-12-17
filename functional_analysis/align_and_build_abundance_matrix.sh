#!/bin/bash

# Remove genes labeled as "hypothetical protein"
awk '/^>/ && toupper($0) ~ /HYPOTHETICAL/ {bool=1}; /^>/ && toupper($0) !~ /HYPOTHETICAL/ {bool=0}; {if (bool==0) print}' all_gene_catalog > gene_catalog_edited

## align raw reads to gene catalog using bowtie2
# build bowtie2 index
bowtie2-build --threads 20 gene_catalog_edited gene_catalog_edited

num_threads=70

f1=${1}
f2=${2}

read1_filename=${f1}
read2_filename=${f2}
name=${read1_filename%_kneaddata*}
fasta_filename=${3} #catalog
index=${4}

bam_filename=${name}
bam_filename=${bam_filename}'.catalog.bam'
echo "Starting Alignment." >&2
bowtie2 \
    -p ${num_threads} \
    -D 20 \
    -R 3 \
    -N 1 \
    -L 20 \
    -i S,1,0.50 \
    -x ${index} \
    --local \
    -q \
    --quiet \
    --mm \
    -1 ${read1_filename} \
    -2 ${read1_filename} \
    | ./samtools-1.9/samtools view -T ${fasta_filename} -b -h -o ${bam_filename} -


# Sort the bam file
echo 'Sorting the bam file' >&2
./samtools-1.9/samtools 'sort' \
    -l 9 \
    -o ${bam_filename%.*}'.sorted.bam' \
    -O bam \
    -@ ${num_threads} \
    ${bam_filename}

# Cleaning up unsorted bam
echo 'Removing unsorted bam' >&2
rm ${bam_filename}
bam_filename=${bam_filename%.*}'.sorted.bam'

# Index the bam
echo 'Indexing the bam file' >&2
bam_index_filename=${bam_filename%.*}'.bai'
./samtools-1.9/samtools 'index' -b ${bam_filename} ${bam_index_filename}

# calculate relative abundances for each gene per sample
# normalizing: divide by total length of genes and number of reads aligned to gene catalog per sample.
for filename in ./*bam; do

	name=${filename%.catalog.sorted.bam}
	echo ${filename}
	echo ${name}
	countsfile=${name}_alignment_data.tsv
	samtools idxstats -@ ${num_threads} ${filename} > ${countsfile}
	head -n -1 ${countsfile} > temp 
	mv temp ${countsfile}

	echo Normalizing "${name}"
	awk '{printf "%.20f\t",$3 / $2}1' ${countsfile} > foo
	mv foo ${countsfile}
	totalaligncount=$(cut -f1 ${countsfile} | awk '{s+=$1} END {print s}')
	awk '{printf "%.20f\t",$1/"'${totalaligncount}'"}1' ${countsfile} > foo
	cut -f1 foo > ${countsfile}

	sed -i 's/0.00000000000000000000/0.0/g' ${countsfile}

done

# concatenate
# to run this, just have a directory with all of your normalized-by-sample data in tsvs sitting there
for filename in *tsv; do
	name=${filename%_alignment_data.tsv} 
	var=$((var+1))
	if [ $var -eq 1 ]; then
		./samtools-1.9/samtools idxstats -@ 16 $name.catalog.sorted.bam | cut -f1 > fornames
		sed -i '1i\	\' fornames 
		cut -f1 fornames > abundance_matrix.csv	
	fi
	echo $name > temp
	cat temp $filename > foo && mv foo ${filename}_temp
	paste abundance_matrix.csv ${filename}_temp > foo && mv foo abundance_matrix.csv
	rm ${filename}_temp 
done
