#!/usr/bin/env python

'''
quickSort.py

This script sorts a fasta file lexically by sequence id.

Usage: quickSort.py in.fa out.fa

    Invoke with -h or --help to print this documentation.

    in.fa   Fasta file containing the sequences to be sorted.
    out.fa  File in which to save the sorted sequences.

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


def main():
    seqs = []
    with open( sys.argv[1], "rU" ) as inFile:
        for s in SeqIO.parse(inFile, "fasta"):
            seqs.append(s)
    with open( sys.argv[2], "w" ) as output:
        SeqIO.write( sorted(seqs,key=lambda r: r.id), output, "fasta" )


if __name__ == "__main__":

    #check if I should print documentation
    q = lambda x: x in sys.argv
    if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
	print __doc__
	sys.exit(0)
        
    #log command line
    logCmdLine(sys.argv)

    main()

