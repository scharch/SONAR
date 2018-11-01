#!/usr/bin/env python3

"""
getFastaFromList.py

This is a simple utility script for efficiently extracting a subset of
      sequences from a large fasta file.

Usage: getFastaFromList.py -f <seqs.fa> [ -o <output.fa> -l <list.txt> ]

Options:
    -f <seqs.fa>     Fasta file containing the sequences to be subsetted.
    -o <output.fa>   Fasta file in which to save extracted sequences. [default: STDOUT]
    -l <list.txt>    Text file containing list of sequence identifiers to extract. [default: STDIN]

Created by Chaim A Schramm on 2015-04-27.
Added streaming from STDIN on 2016-06-10.
Added streaming to STDOUT on 2017-03-19.
Edited to use Py3 and DocOpt by CAS 2018-08-28.
Added FastQ and GZip support by CAS 2018-11-01.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys, fileinput, gzip
from docopt import docopt
from functools import partial
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def loadAndAnnotate(seqFile, saveDict):
	good = 0

	#fasta or fastq
	read_format="fasta"
	if re.search("\.(fq|fastq)", seqFile) is not None:
		read_format = "fastq"

	#gzip?
	if re.search("gz$", seqFile):
		_open = partial(gzip.open,mode='rt')
	else:
		_open=open
		
	with _open(seqFile) as seqs:
		for s in SeqIO.parse(seqs, read_format):
			if s.id not in saveDict:
				continue
			else:
				good += 1
				s.description = s.description + " " + saveDict.get(s.id, "")
				yield s
				if good % 100000 == 0:
					sys.stderr.write("Loaded %d so far...\n" % good)

def main():

    toSave = dict()
    
    for line in fileinput.input(arguments['-l']):
	    # use id as the key and the rest of the line as value to be added to fasta def line
	    fields = line.strip().split()
	    if len(fields) == 0:
		    continue
	    toSave[ fields[0] ] = " ".join( fields[1:] )

    if arguments['-o'] != "STDOUT":
	    sys.stdout = open(arguments['-o'], "w")
    SeqIO.write(loadAndAnnotate(arguments['-f'], toSave), sys.stdout, "fasta")


if __name__ == '__main__':

	arguments = docopt(__doc__)

	if arguments['-l'] == "STDIN":
		arguments['-l'] = "-"
	
	#log command line
	logCmdLine(sys.argv)

	main()

