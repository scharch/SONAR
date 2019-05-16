#!/usr/bin/env python3

"""
1.5-single_cell_statistics.py

This script parses V(D)J annotations in rearrangements.tsv in the context of the `cell_id`
    field to collapse redundancy within each cell and generate useful summary statistics.
    Unfortunately, I do not currently have a good heuristic for identifying contaminating
    transcript that are present in multiple cells without flagging high-frequency
    recombinations (eg VK3-20), especially with low SHM.
Basic statistics are logged in output/logs/cell_processing.log; detailed per-cell
    information gets saved in output/tables/<project>_cells_stats.tsv; and an updated
    rearrangements file with redundancy and irrelevant cells (see the --save option)
    removed is created as <input_rearrangements>_single-cell.tsv.

Usage: 1.5-single_cell_statistics.py [ --rearrangements rearrangments.tsv --save all ]

Options:
    --rearrangements rearrangements.tsv   AIRR-formatted rearrangments file to be processed. 
                                              [default: output/tables/<project>_rearrangements.tsv]
    --save all                            Which cells to include in the reduced rearrangements
                                              table. Options are 'all', 'good' (excludes cells with
                                              multiple productive chains from the same locus), and
                                              'paired' (includes only cells with exactly 1 productive
                                              heavy chain and 1 productive light chain).
                                              [default: all]

Created by Chaim A Schramm on 2019-03-7.
Added output filtering by CA Schramm 2019-04-09.

Copyright (c) 2019 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, re
import subprocess
import airr
from docopt import docopt
from collections import defaultdict, Counter
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *


	
def main():
	cells_raw = defaultdict( dict )
	most_used = defaultdict( list )
	output = open("%s/%s_cell_stats.tsv"%(prj_tree.tables,prj_name), 'w')
	output.write("cell\tstatus\tisotype\tproductive_IGH\ttotal_IGH\tIGH_junctions\tproductive_IGK\ttotal_IGK\tIGK_junctions\tproductive_IGL\ttotal_IGL\tIGL_junctions\n")
	data = airr.read_rearrangement(arguments['--rearrangements'])
	cells_only = airr.derive_rearrangement(re.sub(".tsv", "_single-cell.tsv", arguments['--rearrangements']), arguments['--rearrangements'])

	#assume cells might not be grouped together, so make a first pass
	#    to collect everything
	for r in data:
		if r['status'] in ['good', 'indel', 'stop', 'nonproductive', 'unique']: #skip irrelevant sequences
			if r['locus'] not in cells_raw[r['cell_id']]:
				cells_raw[r['cell_id']][r['locus']] = [ r ]
			else:
				cells_raw[r['cell_id']][r['locus']].append( r )
			#need better heuristic for this, omit for now
			#if r['cell_id'] not in most_used[r['centroid']]:
			#	most_used[r['centroid']].append(r['cell_id'])

	#now go back and process each cell
	status_list  = [ 'canonical_pair', 'possible_inclusion', 'heavy_only', 'light_only', 'multi_light', 'multi_heavy', 'probable_multiplet', 'none_productive' ]
	status_count = dict( zip( status_list, [0,0,0,0,0,0,0,0] )  )
	status_dict  = dict( )

	for c in cells_raw:
		cell_processed	= defaultdict( list )
		cell_productive = defaultdict( list )
		for locus in cells_raw[c]:
			#Start with the one with the most UMIs
			for rep in sorted( [ r for r in cells_raw[c][locus] ], key=lambda k: k['duplicate_count'], reverse=True ):
				#check if this is a duplicate of a previously kept read
				keep = True
				for previous in cell_processed[locus]:
					#shortcut: assume identical junctions means duplicates
					if previous['junction_aa'] == rep['junction_aa']:
						keep = False
						break
					#heuristic (for 10x data as of March 2019):  omit gaps and cut off possible noise at 5' end
					else:
						cov, score = scoreAlign( quickAlign(previous['sequence_alignment'],rep['sequence_alignment']), countInternalGaps=False, skip=50 )
						if score >= 0.95:
							keep = False
							break
					
				if keep:
					cell_processed[locus].append( rep )
					if rep['status'] == "good": cell_productive[locus].append( rep )

		status = ""
		h_type = ""
		if len(cell_productive['IGH']) == 0:
			if len(cell_productive['IGK']) + len(cell_productive['IGL']) == 0:
				status = "none_productive"
			elif len(cell_productive['IGK']) + len(cell_productive['IGL']) == 1:
				status = "light_only"
			elif len(cell_productive['IGK']) + len(cell_productive['IGL']) > 1:
				status = "multi_light"
		elif len(cell_productive['IGH']) == 1:
			h_type = re.sub("\*.+", "", cell_productive['IGH'][0]['c_call'])
			if len(cell_productive['IGK']) + len(cell_productive['IGL']) == 0:
				status = "heavy_only"
			elif len(cell_productive['IGK']) + len(cell_productive['IGL']) == 1:
				status = "canonical_pair"
			elif len(cell_productive['IGK']) + len(cell_productive['IGL']) == 2:
				status = "possible_inclusion"
			elif len(cell_productive['IGK']) + len(cell_productive['IGL']) > 2:
				status = "probable_multiplet"
		elif len(cell_productive['IGH']) > 1:
				status = "multi_heavy"

		status_count[status] += 1
		status_dict[c] = status

		#print to filtered rearrangements file
		if status in arguments['--save']:
			for loc in cell_processed:
				for chain in cell_processed[loc]:
					cells_only.write( chain )

		#now log the cell
		print( "\t".join( [c, status, h_type,
				   str(len(cell_productive['IGH'])), str(len(cell_processed['IGH'])), ";".join([chain['junction_aa'] for chain in cell_processed['IGH']]),
				   str(len(cell_productive['IGK'])), str(len(cell_processed['IGK'])), ";".join([chain['junction_aa'] for chain in cell_processed['IGK']]),
				   str(len(cell_productive['IGL'])), str(len(cell_processed['IGL'])), ";".join([chain['junction_aa'] for chain in cell_processed['IGL']]) ] ),
		       file=output)


	output.close()

	with open("%s/cell_processing.log"%prj_tree.logs, "w") as log:
		print("\t".join(status_list), file=log)
		print("\t".join([str(status_count[s]) for s in status_list]), file=log)
	
	print("\t".join(status_list))
	print("\t".join([str(status_count[s]) for s in status_list]))

#	print("--------------------------------------------------------------------------------------------")
#	print("id\tnum_cells\t" + "\t".join(status_list))
#	for big in sorted( most_used, key=lambda k: len(most_used[k]), reverse=True ):
#		if len(most_used[big]) < 10:
#			break #just want to see the *really* big ones
#		byCategory = Counter([status_dict[c] for c in most_used[big]])
#		print( "%s\t%d\t%s" % (big, len(most_used[big]), "\t".join( [str(byCategory[s]) for s in status_list] )) )
	
	

if __name__ == '__main__':
	
	arguments = docopt(__doc__)
			
	#load saved locus information
	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)
	
	arguments['--rearrangements'] = re.sub("<project>", prj_name, arguments["--rearrangements"])
	
	if not os.path.isfile(arguments['--rearrangements']):
		sys.exit("Cannot find rearrangments file %s" % arguments['--rearrangements'])
	elif not airr.validate_rearrangement(arguments['--rearrangements']):
		sys.exit("File %s is not in valid AIRR format." % arguments['--rearrangements'])

	if re.match( arguments['--save'], 'all', re.I ):
		arguments['--save'] = [ "canonical_pair", "heavy_only", "light_only", "possible_inclusion", "multi_light", "multi_heavy", "probable_multiplet" ]
	elif re.match( arguments['--save'],'good', re.I ):
		arguments['--save'] = [ "canonical_pair", "heavy_only", "light_only", "possible_inclusion" ]
	elif re.match( arguments['--save'],'paired', re.I ):
		arguments['--save'] = [ "canonical_pair" ]
	else:
		sys.exit( "%s is not a valid option for --save. Please select from 'all', 'good', or 'paired' only." % arguments['--save'] )
	
	#log command line
	logCmdLine(sys.argv)

		
	main()

