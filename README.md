# ancient-microbiome-denovo

Scripts used in the paper titled "Reconstruction of ancient microbial genomes from the human gut". Each folder contains its own README.

## Abstract
Loss of gut microbial diversity in modern-day industrial populations is associated with  39 chronic diseases, underscoring the importance of studying our ancestral gut microbiome. In this study,  40 we performed the first large-scale de novo assembly of microbial genomes from ancient microbiome  41 samples. From eight authenticated human paleofeces (1,000-2,000 years old) with exceptionally well 42 preserved DNA from the Southwestern U.S. and Mexico, we reconstructed 498 medium- and high 43 quality genomes. Among the 181 genomes with the strongest evidence of being ancient and of human  44 gut origin, 39% represent novel species-level genome bins (SGBs). Phylogenetic analysis and tip dating  45 suggest an approximate diversification timeline for the key human symbiont Methanobrevibacter  46 smithii. Moreover, compared to 789 modern human gut microbiome samples from eight countries, our  47 data shows that this ancient gut microbiome is more similar to the modern-day non-industrialized human  48 gut microbiome relative to the modern-day industrialized gut microbiome. Functional profiling of the  49 paleofeces reveals significantly lower abundance of antibiotic resistance and mucin-degrading genes and  50 enrichment of mobile genetic elements relative to the modern-day industrial gut microbiome. This work  51 opens the door to discovering novel gut microbes from ancient microbiome samples and interrogating  52 the evolutionary history of the human gut microbiota by genome reconstruction from paleofeces. 

## Directories
### aDNA_damage
Contains scripts to assess ancient DNA (aDNA) damage levels.

### de_novo_assembly
Scripts to de novo assemble raw reads into contigs, bin contigs into MAGs, assess the quality and novelty of the MAGs, dereplicate MAGs, and assign taxonomy to MAGs.

### functional_analysis
Pipeline to predict genes and CAZymes and process gene catalog.

### phylogenetic_analysis
Pipeline for building phylogenetic trees.

### read_processing_and_quality_control
Codes used for processing raw reads and for quality control.

### reference_based_taxonomy
Reference-based taxonomic analysis of the metagenomes.

### statistics_and_figures
Statistical analyses and scripts used to generate figures in the paper.


