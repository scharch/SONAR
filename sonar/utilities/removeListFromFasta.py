#!/usr/bin/env python3

"""
removeListFromFasta.py

This is a simple utility script for efficiently removing a subset of
      sequences from a large fasta file.

Usage: removeListFromFasta.py -l LIST -f SEQS -o OUTPUT

Options:
    -l LIST     Text file containing list of sequence identifiers to extract
    -f FASTA    Fasta file containing the sequences to be subsetted
    -o OUTPUT   Fasta file in which to save extracted sequences

Created by Chaim A Schramm on 2015-04-27.
Edited to use Py3 and DocOpt by CAS 2018-08-29

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
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
					print( "Loaded %d so far..." % good )

def main():
	
	with open(arguments['-l'], "rU") as handle:
		remove = [ line.strip().split()[0] for line in handle if line.strip() != ""]
		
	with open(arguments['-o'], "w") as output:
		SeqIO.write(loadAndAnnotate(arguments['-f'], remove), output, "fasta")


if __name__ == '__main__':

	arguments = docopt(__doc__)
	
	#log command line
	logCmdLine(sys.argv)

	main()

