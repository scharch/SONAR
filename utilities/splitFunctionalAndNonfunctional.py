#!/usr/bin/env python3

"""
splitFunctionalAndNonfunctional.py

This script parses the status of each read in a lineage to separate
   out those that are functional and those that are from passenger alleles.

Usage: splitFunctionalAndNonfunctional.py [ -a lineageNotations.fa -f functionalLineages.fa -n nonfunctional.fa]

Options:
    -a --all lineageNotations.fa          All sequences with lineage notations.
                                             [default: output/sequences/nucleotide/<project>_allJ_unique_lineageNotations.fa]
    -f --func functionalLineages.fa       Where to save the representatives of functional lineages.
                                             [default: output/sequences/nucleotide/<project>_functionalLineages.fa]
    -n --nonf nonfunctional.fa            Where to save the representatives of nonfunctional lineages.
                                             [default: output/sequences/nucleotide/<project>_nonFunctionalLineages.fa]

Created by Larissa K. Ault 2018-07-11.
Renamed, added functionality, and documented by Chaim A. Schramm 2018-08-30.
Altered to run before 5.1 and discard bad reads from functional lineages by CAS 2018-08-31.
Added lineage size correction for functional lineages with discarded bad reads by CAS 2018-09-06.
Updated for AIRR-format compatibility by CAS 2018-10-18.
Made sequence parsing regex more flexible by CAS 2018-11-14.

Copyright (c) 2011-2018 Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

"""

import sys, re
from docopt import docopt
from Bio import SeqIO
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def fixLineageSize( seq ):
	rep = re.search("clone_rep=(\\d+)", seq.description)
	seq.description = re.sub( "clone_count=\d+", "clone_count=%d"%lineages[rep.group(1)]['size'], seq.description )
	return seq


def main():

	global lineages
	
	#parse lineages
	lineages = dict()
	reads	 = dict()
	
	for sequence in SeqIO.parse( open(arguments['--all'], "r"), "fasta" ):

		info = re.search("status=(?P<status>\\S+).*(?:cluster_count=(?P<cluster>\\d+))? clone_id=\d+ clone_rep=(?P<rep>\\d+) clone_count=(?P<count>\\d+)", sequence.description)

		if info:
			
			reads[sequence.id] = { 'status':info.group('status'), 'rep':info.group('rep') }

			if info.group('rep') not in lineages:
				if info.group('status') == "nonproductive":
					lineages[info.group('rep')] = { 'status':"nonproductive", 'size':int(info.group('count')) }
				else:
					lineages[info.group('rep')] = { 'status':"functional", 'size':int(info.group('count')) }

			if info.group('status') == "nonproductive":
				if lineages[info.group('rep')]['status'] == "functional":
					lineages[info.group('rep')]['status'] = "discard"
			else:
				if lineages[info.group('rep')]['status'] == "nonproductive":
					lineages[info.group('rep')]['status'] = "discard"
				elif info.group('status') != "good":
					if info.group('cluster') is not None:
						lineages[info.group('rep')]['size'] -= int(info.group('cluster'))
					else:
						lineages[info.group('rep')]['size'] -= 1
		else:
			#shouldn't happen, but this prevents unexplained KeyErrors from crashing the script
			print("Warning, could not find status and lineage annotations for sequence %s" % sequence.id)
			reads[sequence.id] = { 'status':None, 'rep':sequence.id }
			lineages[sequence.id] = { 'status':"discard" }

		
	#now do output
	with open(arguments['--func'], "w") as fHandle:
		SeqIO.write( [ fixLineageSize(s) for s in SeqIO.parse(open(arguments['--all'], "r"), "fasta") if reads[s.id]['status']=="good" and lineages[reads[s.id]['rep']]['status']=="functional" ], fHandle, "fasta" )

	with open(arguments['--nonf'], "w") as nHandle:
		SeqIO.write( [ s for s in SeqIO.parse(open(arguments['--all'], "r"), "fasta") if lineages[reads[s.id]['rep']]['status']=="nonproductive" ], nHandle, "fasta" )

	

if __name__ == "__main__":

	#parse arguments	
	arguments = docopt(__doc__)

	for arg in arguments:
		arguments[arg] = re.sub("<project>",fullpath2last_folder(os.getcwd()),arguments[arg])

	#log command line
	logCmdLine(sys.argv)

	main()
