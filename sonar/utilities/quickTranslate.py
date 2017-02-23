#!/usr/bin/python

'''
quickTranslate.py nt_in.fa aa_out.fa
'''

import sys
from Bio import SeqIO


def loadAndTranslate():
    with open( sys.argv[1], "rU" ) as inFile:
        for seq in SeqIO.parse(inFile, "fasta"):
            seq.seq = seq.seq.translate()
            yield seq


def main():
    with open( sys.argv[2], "w" ) as output:
        SeqIO.write( loadAndTranslate(), output, "fasta" )

if __name__ == "__main__":
    main()

