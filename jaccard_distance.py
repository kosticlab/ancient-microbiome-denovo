#!/usr/bin/env python3
import numpy as np
import csv
from numpy import genfromtxt
from skbio.diversity.beta import pw_distances
taxonomy_tbl = np.genfromtxt("BAM_9ancient_hmp_protein.csv", delimiter=',', dtype=float, skip_header=1, usecols=range(1,157))
taxonomy_trans = taxonomy_tbl.transpose()

with open (r'BAM_9ancient_hmp_protein.csv') as csvfile:
    csv_reader = csv.reader(csvfile)
    sample_ID = next(csv_reader)
sample_ID.pop(0)
print (sample_ID[1])

#jaccard_distance
j_dm = pw_distances(taxonomy_trans, sample_ID, "jaccard")
print (j_dm[0:3])
j_dm.write('jaccard__9ancient_hmp_protein.csv')
