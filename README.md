# ancient-microbiome-denovo

Scripts used in the paper titled "Reconstruction of ancient microbial genomes from the human gut". Each folder contains its own README.

## Abstract
Loss of gut microbial diversity in industrial populations is associated with chronic diseases, underscoring the importance of studying our ancestral gut microbiome. However, relatively little is known about the composition of pre-industrial gut microbiomes. In this study, we performed the first large-scale de novo assembly of microbial genomes from paleofeces. From eight authenticated human paleofeces (1,000-2,000 years old) with well-preserved DNA from the Southwestern U.S. and Mexico, we reconstructed 498 medium- and high-quality microbial genomes. Among the 181 genomes with the strongest evidence of being ancient and of human gut origin, 39% represent novel species-level genome bins. Tip dating suggests an approximate diversification timeline for the key human symbiont Methanobrevibacter smithii. Compared to 789 present-day human gut microbiome samples from eight countries, the paleofeces are more similar to non-industrialized than industrialized human gut microbiomes. Functional profiling of the paleofeces reveals significantly lower abundance of antibiotic resistance and mucin-degrading genes, as well as enrichment of mobile genetic elements relative to industrial gut microbiomes. This work opens the door to discovering and characterizing novel gut microbes from ancient microbiomes and interrogating the evolutionary history of the human gut microbiota through genome reconstruction from paleofeces. 

## Directories
### aDNA_damage
Contains scripts to assess ancient DNA (aDNA) damage levels.

### de_novo_assembly
Scripts to de novo assemble raw reads into contigs, bin contigs into MAGs, assess the quality and novelty of the MAGs, dereplicate MAGs, and assign taxonomy to MAGs.

### functional_analysis
Pipeline to predict genes and CAZymes and process gene catalog.

### Msmithii_molecular_clocking
Codes used to perform molecular clocking for Methanobrevibacter smithii

### parasite_analysis
Scripts to analyze and visualize parasites identified in the samples.

### phylogenetic_analysis
Pipeline for building phylogenetic trees.

### read_processing_and_quality_control
Codes used for processing raw reads and for quality control.

### reference_based_taxonomy
Reference-based taxonomic analysis of the metagenomes.

### statistics_and_figures
Statistical analyses and scripts used to generate figures in the paper.


