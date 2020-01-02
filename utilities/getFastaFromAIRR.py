#!/usr/bin/env python3

"""
getFastaFromAIRR.py

This is a utility script to extract sequences from an AIRR-formatted TSV and
      save them as a fasta file for use with eg alignment or phylogenetic tools.

Usage: getFastaFromAIRR.py [ --rearrangements AIRR.tsv --column status... (--regex good... | --list ids.txt...) --output sequences.fa --sequence trim --aa]

Options:
    --rearrangements AIRR.tsv    An AIRR-formatted rearrangements file.
                                    [default: output/tables/<project>_rearrangements.tsv]
    --column sequence_id        Which column of the AIRR tsv to filter on. Can be specified
                                    multiple times for multiple filters, in which case the
                                    output will be the subset of sequences passing *all*
                                    filters. [default: status]
    --regex good                A regular expression to be matched by the value of `column`
                                    in each row. Can be specified multiple times for
                                    multiple filters, in which case the output will be the
                                    subset of sequences passing *all* filters. Defaults to
                                    no filter (save all sequences).
    --list ids.txt              Instead of a pattern to match, a list of explicit values
                                    (eg sequence_ids) can be specified in a text file, one
                                    value per line. Multiple files can be specified for
                                    filtering different columns, as for --regex. A single
                                    list can also be specified via STDIN instead of a file
                                    name.
    --output sequences.fa       File in which to save extracted ids or sequences.
                                    [default: STDOUT]
    --sequence trim             Which sequence to save. Options are raw, trim, and junction.
                                    [default: trim]
    --aa                        Flag to request output in amino acids instead of nucleotides.
                                    [default: False]

Created by Chaim A Schramm on 2020-01-02.

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

	annotationList = []
	for ind,col in enumerate(arguments['--column']):
		annotationList.append( {'column':col, 'list':[]} )
		if len(arguments['--list']) > 0:
			for line in fileinput.input( arguments['--list'][ind] ):
				annotationList[ind]['list'].append(line.strip())
		else:
			annotationList[ind]['list'].append(arguments['--regex'][ind])

	if arguments['--output'] != "STDOUT":
		sys.stdout = open(arguments['--output'], "w")

	field = 'sequence_alignment'
	if arguments['--sequence'] == 'raw':
		field = 'sequence'
	elif arguments['--sequence'] == 'junction':
		field = 'junction'

	SeqIO.write(airrToFasta(filterAirrTsv(arguments['--rearrangements'], annotationList),field=field,aa=arguments['--aa']), sys.stdout, "fasta")



if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)
	arguments['--rearrangements'] = re.sub("<project>", prj_name, arguments["--rearrangements"])

	if not os.path.isfile(arguments['--rearrangements']):
		sys.exit(f"Cannot find specified rearrangement file {arguments['--rearrangements']}")

	if not arguments['--sequence'] in ['raw', 'trim', 'junction']:
		sys.exit("Error: options for --sequence are raw, trim, or junction only.")

	if len(arguments['--list'])==0 and len(arguments['--regex'])==0:
		arguments['--regex'].append(".")

	if len(arguments['--list'])==1 and arguments['--list'][0] == "STDIN":
		arguments['--list'][0] = "-"

	#kludge to make sure everything is aligned properly
	if len(arguments['--column']) > 1:
		if max( len(arguments['--regex']),len(arguments['--list']) ) > 1:
			if not len(arguments['--column']) == max( len(arguments['--regex']),len(arguments['--list']) ):
				sys.exit("Unequal number of columns and filters!")
		else:
			#one filter but multiple columns
			if len(arguments['--list']) > 0:
				arguments['--list'] = arguments['--list'] * len(arguments['--column'])
			else:
				arguments['--regex'] = arguments['--regex'] * len(arguments['--column'])
	elif max( len(arguments['--regex']),len(arguments['--list']) ) > 1:
		#multiple filters but one columns
		arguments['--column'] = arguments['--column'] * max( len(arguments['--regex']),len(arguments['--list']) )

	#log command line
	logCmdLine(sys.argv)

	main()
