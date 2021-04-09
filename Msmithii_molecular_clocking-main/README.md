# Msmithii_molecular_clocking
Code available for Marsha et al. 2021 - "Reconstruction of ancient microbial genomes from the human gut"

# Optimizing molecular clocking analysis by inspecting sites affected by ancient DNA damage and strain heterogeneity

### Requirements:
* Python >= 3.6
* biopython 
* pysam >= 0.15.4

## Detecting artificial mutations in the MSA from ancient DNA damage effects

~~~Bash
ancient_artefacts_detector.py artefacts_detecting_example_msa.fna \
--output_CT CT_artefacts.txt \
--output_GA GA_artefacts.txt \ 
~~~

Where: 
`<artefacts_detecting_example_msa.fna>` The initial MSA for molecular clocking analysis
`<CT_artefacts.txt>` Detected artificial sites of C/T substitution.
`<GA_artefacts.txt>` Detected artificial sites of G/A substitution. 

Outputs:

two one-column text files (for CT and GA artificial substitutions respectively) containing detected sites' coordinates in the MSA file. Please check `example_data`

~~~Bash
usage: ancient_artefacts_detector.py [-h] [-opt_CT OUTPUT_CT]
                                     [-opt_GA OUTPUT_GA]
                                     [-a_ratio ANCIENT_RATIO]
                                     [-m_ratio MODERN_RATIO]
                                     [input_msa]

positional arguments:
  input_msa             Input a multiple sequence alignment in fasta form.

optional arguments:
  -h, --help            show this help message and exit
  -opt_CT OUTPUT_CT, --output_CT OUTPUT_CT
                        Specify the output file name for suspicious CT
                        artefacts.
  -opt_GA OUTPUT_GA, --output_GA OUTPUT_GA
                        Specify the output file name for suspicious GA
                        artefacts.
  -a_ratio ANCIENT_RATIO, --ancient_ratio ANCIENT_RATIO
                        The ratio of T or A to other nucleotides among ancient
                        sequences. From 0 to 1.
  -m_ratio MODERN_RATIO, --modern_ratio MODERN_RATIO
                        The ratio of C or G to other nucleotides among modern
                        sequences. From 0 to 1.

~~~

## Detecting sites in the MSA from multiple strains

~~~Bash
strain_heterogeneity_detector.py \
strain_heterogeneity_detecting_example.fna \
strain_heterogeneity_detecting_example.bam \
m__ppa3_SchirmerM_2016__G88704__bin.11 > heterogeneity_sites.txt
~~~

Where:
`<strain_heterogeneity_detecting_example.fna>` A FAST file containing one single continuous nucleotide sequence extracted from MSA
`<strain_heterogeneity_detecting_example.bam>` A BAM file containing reads-contig alignment by mapping metagenomic reads against the extracted sequence from MSA.
`<m__ppa3_SchirmerM_2016__G88704__bin.11>` The header name of the contig in `<strain_heterogeneity_detecting_example.fna>`  

Output:

A four-column text file: First column indicates the detected sites' coordinate in the MSA; Second column gives the dominant base; Third column is the possible alternatives based on mapping analysis; Fourth column is the dominance alllel rate for the detected sites. Note: here just reported sites < 0.8 dominance rate. Please check `example_results`

~~~Bash
usage: strain_heterogeneity_detector.py [-h] [mag] [bam] [contigs]

positional arguments:
  mag         The mag or genome in fasta file.
  bam         reads to mag alignment bam file.
  contigs     Specify the header of reference used in reads-contig mapping.

optional arguments:
  -h, --help  show this help message and exit

~~~

Note: all example data and respective results can be found in folder `example_data` and `example_results`
