## Overview
This directory contains the scripts used to assess ancient DNA damage patterns in the metagenomes. Damage patterns of microbial DNA were assessed using damageprofiler.sh, while damage patterns of human DNA were performed using mapdamage.sh.

## Files
### damageprofiler.sh
This script was used to assess microbial DNA damage patterns in the metagenomes. Each of the medium-quality and high-quality reconstructed genomes was used as reference for its respective sample. For each genome, reads were mapped to each contig, DamageProfiler was run per contig, and the average damage levels and damage variation across reads per contig were calculated.

### mapdamage.sh
This script was used to assess human DNA damage patterns in the metagenomes. Reads were first mapped to the human mtDNA reference genome (rCRS), and the alignment files were used as input for mapDamage2.0.
