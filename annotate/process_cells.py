#!/usr/bin/env python3

"""
process_cells.py

This is a helper script to speed up single cell mode.

Usage: process_cells.py PICKLE

Options:
    PICKLE          Pickled umi dictionary produced by 1.5-single_cell_statistics.py

Split out from 1.5-single_cell_statistics.py by Chaim A Schramm on 2024-06-24.

Copyright (c) 2019-2024 Vaccine Research Center, National Institutes of Health, USA.
    All rights reserved.

"""

import sys, os, re, pickle
from docopt import docopt
import datetime
from collections import defaultdict
from Bio import SeqIO

try:
    from SONAR.annotate import *
except ImportError:
    find_SONAR = sys.argv[0].split("SONAR/annotate")
    sys.path.append(find_SONAR[0])
    from SONAR.annotate import *


def main():

	with open(arguments['PICKLE'], 'rb') as pickle_in:
 		chunk_data = pickle.load(pickle_in)

	cell_processed  = {}
	cell_productive	= {}

	print( f"{datetime.datetime.now()}: processing cells from {arguments['PICKLE']}..." )

	for c, cells_raw in chunk_data:
		cell_processed[c]	= defaultdict( list )
		cell_productive[c] = defaultdict( list )
		for locus in cells_raw:
			#Start with the one with the most UMIs
			for rep in sorted( [ r for r in cells_raw[locus] ], key=lambda k: k['duplicate_count'] or 0, reverse=True ):
				#check if this is a duplicate of a previously kept read
				keep = True
				for previous in cell_processed[c][locus]:
					#shortcut: assume identical junctions means duplicates
					if previous['junction_aa'] == rep['junction_aa']:
						keep = False
						if previous['duplicate_count'] is not None and rep['duplicate_count'] is not None: previous['duplicate_count'] += rep['duplicate_count']
						if previous['consensus_count'] is not None and rep['duplicate_count'] is not None: previous['consensus_count'] += rep['consensus_count']
						break
					#heuristic (for 10x data as of March 2019):  omit gaps and cut off possible noise at 5' end
					else: #elif len(cells_raw[locus])<5: #temporary speed up so I can go to bed.
						score, cov = scoreAlign( quickAlign(previous['sequence_alignment'],rep['sequence_alignment']), countInternalGaps=False, skip=50 )
						if score >= 0.95:
							keep = False
							if previous['duplicate_count'] is not None and rep['duplicate_count'] is not None: previous['duplicate_count'] += rep['duplicate_count']
							if previous['consensus_count'] is not None and rep['duplicate_count'] is not None: previous['consensus_count'] += rep['consensus_count']
							break

				if keep:
					cell_processed[c][locus].append( rep )
					if rep['status'] == "good": cell_productive[c][locus].append( rep )

	with open(re.sub("cells_raw","cells_processed",arguments["PICKLE"]), 'wb') as pickle_out:
		pickle.dump( {'processed':cell_processed, 'productive':cell_productive}, pickle_out )

	print( f"{datetime.datetime.now()}: Finished processing {arguments['PICKLE']}" )


if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree = ProjectFolders( os.getcwd() )
	prj_name = fullpath2last_folder(prj_tree.home)

	#log command line
	#logCmdLine(sys.argv)


	main()
