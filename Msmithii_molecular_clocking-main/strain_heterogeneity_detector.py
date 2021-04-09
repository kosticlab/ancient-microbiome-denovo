#!/usr/bin/env python

import sys
import argparse
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

"""
NAME: strain_heterogeneity_detector.py

DESCRIPTION: strain_heterogeneity_detector.py is a python script which takes reads-contigs bam
			 file inputs. It screens all contigs site-by-site, and spots out those stemming from
			 multiple strains (allele dominance rate < 0.8). 
"""

__author__ = "Kun D. Huang"
__version__ = "0.1"
__date__ = "18.06.2020"



def read_args(args):

	parser = argparse.ArgumentParser()

	parser.add_argument('m',
					nargs = "?",
					metavar = "mag",
					help = 'The mag or genome in fasta file.',
					type = str)
	
	parser.add_argument('b',
						nargs = "?",
						metavar = "bam",
						help = 'reads to mag alignment bam file.',
						type = str)
	parser.add_argument('c',
						nargs = "?",
						metavar = "contigs",
						help = 'Specify the header of reference used in reads-contig mapping.',
						type = str)	


	return vars(parser.parse_args())



class bam_miner(object):
	"""
	object bam_miner is to mining bam file.
	"""
	def __init__(self, contig, bam):
		self.bam = bam
		self.contig = contig

	def mpileup(self):
		bamfile = pysam.AlignmentFile(self.bam, "rb" )
		pos_sub_dict = {}
		for pileupcolumn in bamfile.pileup(self.contig):
			# print("\ncoverage at base %s = %s" %
			# 	(pileupcolumn.pos, pileupcolumn.n))
			altering_dict = {"A":0, "T":0, "G":0, "C":0}
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					# query position is None if is_del or is_refskip is set.
					pos = pileupcolumn.pos
					read_num = pileupcolumn.n
					nucleotide = pileupread.alignment.query_sequence[pileupread.query_position]
					# print(pos, read_num, nucleotide)
					if nucleotide not in altering_dict:
						altering_dict[nucleotide] = 1
					else:
						altering_dict[nucleotide] += 1
			pos_sub_dict[pileupcolumn.pos] = altering_dict
					# print ('\tbase in read %s = %s at %s' %
					# (pileupread.alignment.query_name,
					# 	pileupread.alignment.query_sequence[pileupread.query_position], pileupread.query_position))

		bamfile.close()
		return pos_sub_dict



if __name__== '__main__':

	pars = read_args(sys.argv)
	obj = bam_miner(pars['c'], pars['b'])
	pos_sus_db = obj.mpileup()
	ref_seq = list(SeqIO.parse(open(pars['m']), "fasta"))[0].seq.upper()
	
	count = 0
	for i in range(1, len(ref_seq)):
		if ref_seq[i] != '-' and i in pos_sus_db:
			T_frq = pos_sus_db[i]['T']
			A_frq = pos_sus_db[i]['A']
			C_frq = pos_sus_db[i]['C']
			G_frq = pos_sus_db[i]['G']
			sum_sub  = T_frq + A_frq + C_frq + G_frq
			if sum_sub != 0:
				dominance = pos_sus_db[i][ref_seq[i]]/sum_sub
				if dominance < 0.8:
					# count += 1
					# print(i+1)
					print(i+1, ref_seq[i], pos_sus_db[i], dominance)
			else:
				print(i+1, ref_seq[i], pos_sus_db[i], 0)
				# print(i+1)




