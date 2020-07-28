#!/usr/bin/env python3

"""
getFastaFromAIRR.py

This is a utility script to extract sequences from an AIRR-formatted TSV and
      save them as a fasta file for use with eg alignment or phylogenetic tools.
      To only save certain sequences, filter the AIRR TSV using utilities/filterAIRR.py
      and pipe the output to this script with `--rearrangements STDIN`.

Usage: getFastaFromAIRR.py [ --rearrangements AIRR.tsv --output sequences.fa --sequence trim --aa]

Options:
    --rearrangements AIRR.tsv    An AIRR-formatted rearrangements file. Use 'STDIN' to get a
                                    stream. [default: output/tables/<project>_rearrangements.tsv]
    --output sequences.fa       File in which to save extracted ids or sequences.
                                    [default: STDOUT]
    --sequence trim             Which sequence to save. Options are raw, trim, and junction.
                                    [default: trim]
    --aa                        Flag to request output in amino acids instead of nucleotides.
                                    [default: False]

Created by Chaim A Schramm on 2020-01-02.
Added --equal option by CA Schramm on 2020-02-21.
Removed filtering options (use filterAIRR.py)  and added streaming
                         input by CA Schramm 2020-07-02.

Copyright (c) 2020 Vaccine Research Center, National Institutes of
                         Health, USA. All rights reserved.

"""

import sys, re, fileinput
from docopt import docopt
from collections import defaultdict
import airr

try:
	from SONAR import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/utilities")
	sys.path.append(find_SONAR[0])
	from SONAR import *



def main():

	if arguments['--output'] != "STDOUT":
		sys.stdout = open(arguments['--output'], "w")

	field = 'sequence_alignment'
	if arguments['--sequence'] == 'raw':
		field = 'sequence'
	elif arguments['--sequence'] == 'junction':
		field = 'junction'


	infile = arguments['--rearrangements']
	if arguments['--rearrangements'] == "STDIN":
		infile = "-"
	reader = airr.io.RearrangementReader(fileinput.input(infile))

	SeqIO.write(airrToFasta(reader, field=field, aa=arguments['--aa']), sys.stdout, "fasta")



if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)
	arguments['--rearrangements'] = re.sub("<project>", prj_name, arguments["--rearrangements"])

	if not os.path.isfile(arguments['--rearrangements']) and arguments['--rearrangements'] != "STDIN":
		sys.exit(f"Cannot find specified rearrangement file {arguments['--rearrangements']}")

	if not arguments['--sequence'] in ['raw', 'trim', 'junction']:
		sys.exit("Error: options for --sequence are raw, trim, or junction only.")

	#log command line
	logCmdLine(sys.argv)

	main()
