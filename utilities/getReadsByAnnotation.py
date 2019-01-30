#!/usr/bin/env python3

"""
getReadsByAnnotations.py

This is a utility script to find all sequences in a large fasta that have
      an annotation of interest (eg a particular V gene assignment or
      CDR3 motif).

Usage: getReadsByAnnotations.py -f file.fa ( -a annotation | -l list ) [ -t seq -o output.fa -m 0 ]

Options:
    -f file.fa      Fasta file containing the sequences to be subsetted.
    -a annotation   Regular expression to look for in the fasta def line
                       of each sequence.
    -l list         A list of annotations to look for. STDIN or a file name.
    -t seq          Type of output. Options are 'seq' (sequences in fasta
                       format) or 'id' (list of sequence ids). [default: seq]
    -o output.fa    File in which to save extracted ids or sequences. 
                       [default: STDOUT]
    -m 0            Max number of matches to find (0 for all matches).
                       [default: 0]


Created by Chaim A Schramm on 2018-11-01.
Added max matches and list options by CAS 2018-11-13.
Added option to output ids only by CAS 2018-12-12.

Copyright (c) 2011-2018 Vaccine Research Center, National Institutes of
                         Health, USA. All rights reserved.

"""

import sys, re, fileinput
from docopt import docopt
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def annotationInList(annotation,refList):
	for a in refList:
		if re.search( a, annotation ):
			return True
	return False

	
def checkAnnotation(seqFile, annotationList):
	good = 0
	read_format="fasta"
	if re.search("\.(fq|fastq)$", seqFile) is not None:
		read_format="fastq"
	with open(seqFile, "r") as seqs:
		for s in SeqIO.parse(seqs, read_format):
			if annotationInList( s.description, annotationList ):
				good += 1
				yield s
				if good == arguments['-m']:
					break
				if good % 100000 == 0:
					sys.stderr.write("Loaded %d so far...\n" % good)



def main():

	annotationList = []
	if arguments['-l'] is not None:
		for line in fileinput.input(arguments['-l']):
			annotationList.append(line.strip())
	else:
		annotationList.append(arguments['-a'])
		
	if arguments['-o'] != "STDOUT":
		sys.stdout = open(arguments['-o'], "w")

	if arguments['-t'] == 'seq':
		SeqIO.write(checkAnnotation(arguments['-f'], annotationList), sys.stdout, "fasta")
	else:
		for seq in checkAnnotation(arguments['-f'], annotationList):
			sys.stdout.write("%s\n"%seq.id)


    
if __name__ == '__main__':

	arguments = docopt(__doc__)
	arguments['-m'] = int(arguments['-m'])

	if arguments['-l'] == "STDIN":
		arguments['-l'] = "-"

	if arguments['-t'] not in ['seq', 'id']:
		sys.exit("Valid values for -t are 'seq' and 'id' only.\n")
	
	#log command line
	logCmdLine(sys.argv)

	main()

