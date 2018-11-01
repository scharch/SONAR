#!/usr/bin/env python3

"""
getReadsByAnnotations.py

This is a utility script to find all sequences in a large fasta that have
      an annotation of interest (eg a particular V gene assignment or
      CDR3 motif).

Usage: getReadsByAnnotations.py -f file.fa -a annotation [ -o output.fa ]

Options:
    -f file.fa      Fasta file containing the sequences to be subsetted.
    -a annotation   Regular expression to look for in the fasta def line
                       of each sequence.
    -o output.fa    Fasta file in which to save extracted sequences. 
                       [default: STDOUT]


Created by Chaim A Schramm on 2018-11-01.

Copyright (c) 2011-2018 Vaccine Research Center, National Institutes of
                         Health, USA. All rights reserved.

"""

import sys, re
from docopt import docopt
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def checkAnnotation(seqFile):
	good = 0
        read_format="fasta"
        if re.search("\.(fq|fastq)$", seqFile) is not None:
                read_format="fastq"
	with open(seqFile, "rU") as seqs:
		for s in SeqIO.parse(seqs, read_format):
			if re.search( arguments['-a'], s.description ):
				good += 1
				yield s
				if good % 100000 == 0:
					sys.stderr.write("Loaded %d so far...\n" % good)



def main():

    if arguments['-o'] != "STDOUT":
	    sys.stdout = open(arguments['-o'], "w")
    SeqIO.write(checkAnnotation(arguments['-f']), sys.stdout, "fasta")



    
if __name__ == '__main__':

	arguments = docopt(__doc__)

	#log command line
	logCmdLine(sys.argv)

	main()

