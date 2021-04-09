#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys
import argparse


"""
NAME: ancient_artegacts_detector.py
DESCRIPTION: metaclock_combiner.py is a python script which inspects artificial mutation in an
             alignment due to post-mortem DNA damage effect.  
"""

__author__ = "Kun D. Huang"
__version__ = "0.1"
__date__ = "21.07.2020"



def read_args(args):

    parser = argparse.ArgumentParser()

    parser.add_argument('ipt_msa',
                        nargs = "?",
                        metavar = "input_msa",
                        help = "Input a multiple sequence alignment in fasta form.",
                        type = str)

    parser.add_argument('-opt_CT',
                        '--output_CT',
                        help = 'Specify the output file name for suspicious CT artefacts.',
                        type = str)

    parser.add_argument('-opt_GA',
                        '--output_GA',
                        help = 'Specify the output file name for suspicious GA artefacts.',
                        type = str)

    parser.add_argument('-a_ratio',
                        '--ancient_ratio',
                        help = 'The ratio of T or A to other nucleotides among ancient sequences. From 0 to 1.',
                        type = float,
                        default = 1.0)

    parser.add_argument('-m_ratio',
                        '--modern_ratio',
                        help = 'The ratio of C or G to other nucleotides among modern sequences. From 0 to 1.',
                        type = float,
                        default = 1.0)


    return vars(parser.parse_args())



def aln_partition(aln):
    """
    This function partitions whole alignment into
    modern AlignIO object and ancient AlignIO obeject  
    """
    m_seq_records, a_seq_records = [], []
    for r in aln:
        if r.name.startswith('a__'):
            a_seq_records.append(SeqRecord(Seq(str(r.seq).upper(), generic_dna), r.name, description = ''))
        elif r.name.startswith('m__'):
            m_seq_records.append(SeqRecord(Seq(str(r.seq).upper(), generic_dna), r.name, description = ''))
        else:
            sys.exit('Please label sequences correctly!')
    
    return MultipleSeqAlignment(a_seq_records), MultipleSeqAlignment(m_seq_records)

def correct_empty_lst(lst):
    if len(lst) > 0:
        return len(lst)
    else:
        return len(range(10000000000)) 

def suspicious_pos_detector(a_aln, m_aln, _type, opt_file, a_ratio = 1.0, m_ratio = 1.0):
    
    m_base = _type[0]
    a_base = _type[1]

    opt = open(opt_file, 'w')
    for i in range(len(a_aln[1,:])):
        col_a_group_nogap = [i for i in a_aln[:,i]]
        col_m_group_nogap = [i for i in m_aln[:,i] if i != '-']
        a_ratio_est = col_a_group_nogap.count(a_base)/len(col_a_group_nogap)
        m_ratio_est = col_m_group_nogap.count(m_base)/correct_empty_lst(col_m_group_nogap)
        if (a_ratio_est >= a_ratio) and (m_ratio_est >= m_ratio):
            opt.write(str(i+1)+'\n')
    opt.close()



if __name__ == "__main__":

    pars = read_args(sys.argv)
    aln_ipt = AlignIO.read(pars['ipt_msa'], 'fasta')
    a_aln_obj = aln_partition(aln_ipt)[0]
    m_aln_obj = aln_partition(aln_ipt)[1]
    
    suspicious_pos_detector(a_aln_obj, m_aln_obj, 'CT' , pars['output_CT'], pars['ancient_ratio'], pars['modern_ratio'])
    suspicious_pos_detector(a_aln_obj, m_aln_obj, 'GA', pars['output_GA'], pars['ancient_ratio'], pars['modern_ratio'])


