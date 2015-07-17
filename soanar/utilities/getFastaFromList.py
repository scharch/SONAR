#!/usr/bin/python

"""
getFastaFromList.py

This is a simple utility script for efficiently extracting a subset of
      sequences from a large fasta file.

Usage: getFastaFromList.py -l list.txt -f seqs.fa -o output.fa

    Invoke with -h or --help to print this documentation.

    l           Text file containing list of sequence identifiers to extract
    f           Fasta file containing the sequences to be subsetted
    o           Fasta file in which to save extracted sequences

Created by Chaim A Schramm on 2015-04-27.
Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

from soanar import *

def main():

    global listFile, inFile, outFile
    
    with open(listFile, "rU") as handle:
        toSave = [ line.strip() for line in handle ]

    seqs = load_seqs_in_dict(inFile, set(toSave))

    with open(outFile, "w") as output:
        SeqIO.write(seqs.values(), output, "fasta")


if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, l="listFile", f="inFile", o="outFile")
        listFile, inFile, outFile = getParas(dict_args, "listFile", "inFile", "outFile")

	main()

