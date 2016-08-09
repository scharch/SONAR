#!/usr/bin/env python

"""
removeListFromFasta.py

This is a simple utility script for efficiently removing a subset of
      sequences from a large fasta file.

Usage: removeListFromFasta.py -l list.txt -f seqs.fa -o output.fa

    Invoke with -h or --help to print this documentation.

    l           Text file containing list of sequence identifiers to extract
    f           Fasta file containing the sequences to be subsetted
    o           Fasta file in which to save extracted sequences

Created by Chaim A Schramm on 2015-04-27.
Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def loadAndAnnotate(seqFile, killList):
	good = 0
	with open(seqFile, "rU") as seqs:
		for s in SeqIO.parse(seqs, "fasta"):
			if s.id in killList:
				continue
			else:
				good += 1
				yield s
				if good % 100000 == 0:
					print "Loaded %d so far..." % good

def main():

    global listFile, inFile, outFile
    
    with open(listFile, "rU") as handle:
	    remove = [ line.strip().split()[0] for line in handle if line.strip() != ""]

    with open(outFile, "w") as output:
        SeqIO.write(loadAndAnnotate(inFile, remove), output, "fasta")


if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	#log command line
	logCmdLine(sys.argv)

	# get parameters from input
	dict_args = processParas(sys.argv, l="listFile", f="inFile", o="outFile")
        listFile, inFile, outFile = getParas(dict_args, "listFile", "inFile", "outFile")

	main()

