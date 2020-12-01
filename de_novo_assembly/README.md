## Overview
This folder contains scripts we used to de novo assemble raw reads into contigs, bin contigs into MAGs, assess the quality and novelty of the MAGs, dereplicate MAGs, and assign taxonomy to MAGs. 

## Files
### assembly_pipeline.sh
De novo assemble raw reads into contigs and bin contigs into MAGs.

### bbmap.sh
Calculate assembly statistics (number of contigs, number of base pairs in contigs, contig N50, contig L50, and the longest contig). 

### checkm.sh
Assess quality (completeness, contamination, genome size in bp, number of contigs, contig N50 values, mean contig length, length of the longest contig, coverage) of the MAGs.

### drep.sh
Calculate pairwise Average Nucleotide Identity (ANI) for the MAGs and dereplicate MAGs of the same species into species-level genome bins (SGBs).

### gtdbtk.sh
Assign taxonomies to the MAGs.

### novelty_assessment.sh
Assess whether each SGB belongs to a known microbial species by calculating Mash distances and Average Nucleotide Identities (ANIs) between the SGBs and the reference genomes. 
