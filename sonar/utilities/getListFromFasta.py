#!/usr/bin/env python3

"""
getListFromFasta.py

This is a simple utility script for getting a list of the sequence identifiers
      from a large fasta file.

Usage: getListFromFasta.py -f <seqs.fa> [ -o <output.list> ]

Options:
    -f <seqs.fa>       Fasta file containing the sequences.
    -o <output.list>   Text file in which to save the list of sequences.
                          [default: STDOUT]

Created by Chaim A Schramm on 2015-04-27.
Added printing to STDOUT on 2016-06-10.
Edited to use Py3 and DocOpt by CAS 2018-08-28.

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


def main():

    if arguments['-o'] != "STDOUT":
	    sys.stdout = open(arguments['-o'], "w")
    for seq in generate_read_fasta(arguments['-f']):
	    print( "%s"%seq.id )


if __name__ == '__main__':

	arguments = docopt(__doc__)

	#log command line
	logCmdLine(sys.argv)

	main()

