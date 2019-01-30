#!/usr/bin/env python3

'''
quickTranslate.py

This script translates a fasta file.

Usage: quickTranslate.py <in.fa> <out.fa>

Options:
    <in.fa>    Fasta file containing the nucleotide sequences to be translated.
    <out.fa>   File in which to save the translated amino acid sequences.

Added to SONAR by Chaim A Schramm 2017-02-24.
Edited to use Py3 and DocOpt by CAS 2018-08-29.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

'''

import sys
from docopt import docopt
from Bio import SeqIO

try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def loadAndTranslate():
	with open( sys.argv[1], "r" ) as inFile:
		for seq in SeqIO.parse(inFile, "fasta"):
			seq.seq = seq.seq.translate()
			yield seq


def main():
	with open( sys.argv[2], "w" ) as output:
		SeqIO.write( loadAndTranslate(), output, "fasta" )

	
if __name__ == "__main__":

	arguments = docopt(__doc__)
	
	#log command line
	logCmdLine(sys.argv)

	main()

