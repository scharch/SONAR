#!/usr/bin/python

'''
quickSort.py nt_in.fa aa_out.fa
'''

import sys
from Bio import SeqIO


def main():
    seqs = []
    with open( sys.argv[1], "rU" ) as inFile:
        for s in SeqIO.parse(inFile, "fasta"):
            seqs.append(s)
    with open( sys.argv[2], "w" ) as output:
        SeqIO.write( sorted(seqs,key=lambda r: r.id), output, "fasta" )


if __name__ == "__main__":
    main()

