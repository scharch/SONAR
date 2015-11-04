#!/usr/bin/env python

"""
getListFromFasta.py

This is a simple utility script for getting a list of the sequence identifiers
      from a large fasta file.

Usage: getListFromFasta.py -f seqs.fa -o output.list

    Invoke with -h or --help to print this documentation.

    f           Fasta file containing the sequences
    o           Text file in which to save the list of sequences

Created by Chaim A Schramm on 2015-04-27.
Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

from sonar import *

def main():

    global inFile, outFile
    with open(outFile, "w") as output:
        [output.write("%s\n"%seq.id) for seq in generate_read_fasta(inFile)]


if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, f="inFile", o="outFile")
        try:
            inFile, outFile = getParas(dict_args, "inFile", "outFile")
        except KeyError:
            print __doc__
            sys.exit(1)


	main()

