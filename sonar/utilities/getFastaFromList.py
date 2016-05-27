#!/usr/bin/env python

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


def loadAndAnnotate(seqFile, saveDict):
	good = 0
	with open(seqFile, "rU") as seqs:
		for s in SeqIO.parse(seqs, "fasta"):
			if s.id not in saveDict:
				continue
			else:
				good += 1
				s.description = s.description + " " + saveDict.get(s.id, "")
				yield s
				if good % 100000 == 0:
					print "Loaded %d so far..." % good

def main():

    global listFile, inFile, outFile
    
    with open(listFile, "rU") as handle:
        #create a dictionary with the id as the key and the rest of the line as value
        #  to be added to fasta def line
	toSave = { k: v for k,v in [(line.strip().split()[0], " ".join(line.strip().split()[1:])) for line in handle if line.strip() != ""] }

    with open(outFile, "w") as output:
        SeqIO.write(loadAndAnnotate(inFile, toSave), output, "fasta")


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

