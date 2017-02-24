#!/usr/bin/env python

'''
quickTranslate.py

This script translates a fasta file.

Usage: quickTranslate.py in.fa out.fa

    Invoke with -h or --help to print this documentation.

    in.fa   Fasta file containing the nucleotide sequences to be translated.
    out.fa  File in which to save the translated amino acid sequences.

Added to SONAR by Chaim A Schramm 2017-02-24.
Copyright (c) 2011-2017 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

'''

import sys
from Bio import SeqIO

try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *



def loadAndTranslate():
    with open( sys.argv[1], "rU" ) as inFile:
        for seq in SeqIO.parse(inFile, "fasta"):
            seq.seq = seq.seq.translate()
            yield seq


def main():
    with open( sys.argv[2], "w" ) as output:
        SeqIO.write( loadAndTranslate(), output, "fasta" )

        
if __name__ == "__main__":

    #check if I should print documentation
    q = lambda x: x in sys.argv
    if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
	print __doc__
	sys.exit(0)
        
    #log command line
    logCmdLine(sys.argv)

    main()

