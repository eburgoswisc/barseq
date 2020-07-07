#!/usr/bin/env python3
"""

Constants that might be used by Barseq

"""

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'


# Linker Sequences (DEFAULT)
LEFT_LINKER = 'ACGAGACGAGCTTCTTATATATGCTTCGCCAG'
RIGHT_LINKER = 'GACTTGACCTGGATGTCTCTACCCACAAGATCG'

# Reverse of linker
LEFT_LINKER_R = 'CTGGCGAAGCATATATAAGAAGCTCGTCTCGT'
RIGHT_LINKER_R = 'CGATCTTGTGGGTAGAGACATCCAGGTCAAGTC'

# Spacers
spacer = 'GCTCATGCACTTGATTCC'

# Illumina Sequencing (For demultiplexing if needed)
adapter_for = 'aatgatacggcgaccaccgagatctacac'
adapter_rev = 'atctcgtatgccgtcttctgcttg'
pad_for = 'tatggtaatt'
pad_rev = 'ctgactgact'