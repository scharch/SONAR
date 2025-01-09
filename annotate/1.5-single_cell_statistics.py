#!/usr/bin/env python3

"""
1.5-single_cell_statistics.py

This script parses V(D)J annotations in rearrangements.tsv in the context of the `cell_id`
    field to collapse redundancy within each cell and generate useful summary statistics.
    It will also pull in cell hashing and feature barcoding information decoded by
    `1.0-preprocess.py` if available.
Unfortunately, I do not currently have a good heuristic for identifying contaminating
    transcripts that are present in multiple cells without flagging high-frequency
    recombinations (eg VK3-20), especially with low SHM.
Basic statistics are logged in output/logs/cell_processing.log; detailed per-cell
    information gets saved in output/tables/<project>_cells_stats.tsv; and an updated
    rearrangements file with redundancy and irrelevant cells (see the --save option)
    removed is created as <input_rearrangements>_single-cell.tsv.

Usage: 1.5-single_cell_statistics.py [ --rearrangements rearrangments.tsv --save all --cluster --threads 1 ]

Options:
    --rearrangements rearrangements.tsv   AIRR-formatted rearrangments file to be processed.
                                              [default: output/tables/<project>_rearrangements.tsv]
    --save all                            Which cells to include in the reduced rearrangements
                                              table. Options are 'all', 'good' (excludes cells with
                                              multiple productive chains from the same locus), and
                                              'paired' (includes only cells with exactly 1 productive
                                              heavy chain and 1 productive light chain).
                                              [default: all]
    --cluster                             Flag to submit chunk jobs to cluster instead of running them
                                              locally. [default: False]
    --threads 1                           Number of threads to use if not running on a cluster.
                                              [default: 1]

Created by Chaim A Schramm on 2019-03-7.
Added output filtering by CA Schramm 2019-04-09.
Added cell hashing and feature barcoding by CA Schramm 2019-10-16.
Added cell_status to output and included read/UMI counts from discarded
                 duplicates by CAS 2019-12-26.
Included nonproductive rearrangments in detection of multiplets by CAS 2020-01-03.

Copyright (c) 2019-2020 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, re, csv
import subprocess
import airr
from docopt import docopt
from collections import defaultdict, Counter
from io import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial
import itertools
import pickle

try:
	from SONAR.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/annotate")
	sys.path.append(find_SONAR[0])
	from SONAR.annotate import *


# a utility function to get us a slice of an iterator, as an iterator
# when working with iterators maximum lazyness is preferred
# from https://stackoverflow.com/a/44502827
def iterator_slice(iterator, length):
	iterator = iter(iterator)
	while True:
		res = tuple(itertools.islice(iterator, length))
		if not res:
			break
		yield res


def processCells(num):
	cmd = f"{SCRIPT_FOLDER}/annotate/process_cells.py {prj_tree.internal}/cells_raw_{num:04}.pickle"
	os.system( cmd )


def main():

	#look for cell hashing
	hashDict = dict()
	sampleList = []
	if os.path.exists(f"{prj_tree.tables}/{prj_name}_hashes.tsv"):
		with open(f"{prj_tree.tables}/{prj_name}_hashes.tsv", 'r') as handle:
			reader = csv.reader(handle, delimiter="\t")
			for row in reader:
				hashDict[row[0]] = [ row[1] ]
				if row[1]=="unknown" or row[1]=="ambiguous":
					continue
				elif not row[1] in sampleList:
					sampleList.append(row[1])
	sampleList.sort()
	sampleList += ["unknown", "ambiguous"]

	#look for feature barcoding
	featureDict = dict()
	if os.path.exists(f"{prj_tree.tables}/{prj_name}_features.tsv"):
		with open(f"{prj_tree.tables}/{prj_name}_features.tsv", 'r') as handle:
			reader = csv.reader(handle, delimiter="\t")
			header = next(reader)
			featureDict["keys"] = header[1:]
			for row in reader:
				featureDict[row[0]] = row[1:]

	cells_raw = defaultdict( dict )
	most_used = defaultdict( list )

	output = open("%s/%s_cell_stats.tsv"%(prj_tree.tables,prj_name), 'w')
	outwriter = csv.writer( output, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE )
	outheader = ["cell","status","isotype"]
	if len(hashDict) > 0:
		outheader += ["hash_sample"]
	if len(featureDict) > 0:
		outheader += featureDict["keys"]
	outheader += ["productive_IGH","total_IGH","IGH_junctions","productive_IGK","total_IGK","IGK_junctions","productive_IGL","total_IGL","IGL_junctions"]
	outwriter.writerow(outheader)

	data = airr.read_rearrangement(arguments['--rearrangements'])
	fields = ["cell_status"]
	if len(hashDict) > 0:
		fields += ["hash_sample"]
	cells_only = airr.derive_rearrangement(re.sub(".tsv", "_single-cell.tsv", arguments['--rearrangements']), arguments['--rearrangements'],fields=fields)

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
	status_count = dict( )
	for sample in sampleList:
		status_count[sample] = dict( zip( status_list, [0,0,0,0,0,0,0,0] )  )
	status_dict  = dict( )

	#generate pickles to pass to consensus algorithm
	dInd = 0
	for chunk in iterator_slice(cells_raw.items(), 10000):
		dInd += 1
		with open( f"{prj_tree.internal}/cells_raw_{dInd:04}.pickle", 'wb') as pickle_out:
			pickle.dump( chunk, pickle_out )

	#spawn subprocesses
	if arguments['--cluster']:
		with open("%s/processCells.sh"%prj_tree.internal, 'w') as jobHandle:
			jobHandle.write(f"#!/bin/bash\n#$ -N processCells\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/process_cells.py {prj_tree.internal}/cells_raw_$NUM.pickle\n\n")
		subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%dInd, "%s/processCells.sh"%prj_tree.internal])
	else:
		partial_cons = partial( processCells )

		pool = Pool(arguments['--threads'])
		blob = pool.map( partial_cons, range(1,dInd+1) )
		pool.close()
		pool.join()

	#collect output
	for p in glob.glob(f"{prj_tree.internal}/cells_processed_*.pickle"):
		with open(p, 'rb') as pickle_in:
			chunk_dict = pickle.load(pickle_in)

			cell_processed  = chunk_dict[ 'processed' ]
			cell_productive = chunk_dict[ 'productive' ]

			for c in cell_processed:
				status = ""
				h_type = ""
				if len(cell_processed[c]['IGH']) > 2 or len(cell_processed[c]['IGK']) > 2 or len(cell_processed[c]['IGL']) > 2:
					status = "probable_multiplet"
				elif len(cell_productive[c]['IGH']) == 0:
					if len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) == 0:
						status = "none_productive"
					elif len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) == 1:
						status = "light_only"
					elif len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) > 1:
						status = "multi_light"
				elif len(cell_productive[c]['IGH']) == 1:
					h_type = re.sub("\*.+", "", cell_productive[c]['IGH'][0]['c_call'])
					if len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) == 0:
						status = "heavy_only"
					elif len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) == 1:
						status = "canonical_pair"
					elif len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) == 2:
						status = "possible_inclusion"
					elif len(cell_productive[c]['IGK']) + len(cell_productive[c]['IGL']) > 2:
						status = "probable_multiplet"
				elif len(cell_productive[c]['IGH']) > 1:
						status = "multi_heavy"

				status_count[ hashDict.get(c,["unknown"])[0] ][ status ] += 1
				status_dict[c] = status


				#print to filtered rearrangements file
				#leave out cells with ambiguous hashing assignments if we are doing any filtering
				if status in arguments['--save']:
					if hashDict.get(c,["unknown"])[0] != "ambiguous" or "probable_multiplet" in arguments['--save']:
						for loc in cell_processed[c]:
							for chain in cell_processed[c][loc]:
								chain['cell_status'] = status
								if len(hashDict)>0:
									chain['hash_sample']=hashDict.get(chain['cell_id'],['unknown'])[0]
								cells_only.write( chain )

				#now log the cell
				outwriter.writerow( [c, status, h_type] + hashDict.get(c,['unknown']*(len(hashDict)>0)) + featureDict.get(c, ['0']*len(featureDict.get("keys",[]))) +
						   [ len(cell_productive[c]['IGH']), len(cell_processed[c]['IGH']), ";".join([chain['junction_aa'] for chain in cell_processed[c]['IGH']]),
						   len(cell_productive[c]['IGK']), len(cell_processed[c]['IGK']), ";".join([chain['junction_aa'] for chain in cell_processed[c]['IGK']]),
						   len(cell_productive[c]['IGL']), len(cell_processed[c]['IGL']), ";".join([chain['junction_aa'] for chain in cell_processed[c]['IGL']]) ] )

	output.close()

	with open("%s/cell_processing.log"%prj_tree.logs, "w") as log:
		print("sample\t" + "\t".join(status_list), file=log)
		print("sample\t" + "\t".join(status_list))
		for sample in sampleList:
			if sum([status_count[sample][s] for s in status_list]) == 0:
				continue #leave out `ambiguous` if it's not a hashed sample
			print("\t".join([sample]+[str(status_count[sample][s]) for s in status_list]), file=log)
			print("\t".join([sample]+[str(status_count[sample][s]) for s in status_list]))


#	print("--------------------------------------------------------------------------------------------")
#	print("id\tnum_cells\t" + "\t".join(status_list))
#	for big in sorted( most_used, key=lambda k: len(most_used[k]), reverse=True ):
#		if len(most_used[big]) < 10:
#			break #just want to see the *really* big ones
#		byCategory = Counter([status_dict[c] for c in most_used[big]])
#		print( "%s\t%d\t%s" % (big, len(most_used[big]), "\t".join( [str(byCategory[s]) for s in status_list] )) )

	#clean up!!
	oldFiles = glob.glob(f"{prj_tree.internal}/cells_*.pickle")
	if len(oldFiles) > 0:
		[os.remove(f) for f in oldFiles]


if __name__ == '__main__':

	arguments = docopt(__doc__)

	#load saved locus information
	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)

	arguments['--rearrangements'] = re.sub("<project>", prj_name, arguments["--rearrangements"])
	arguments['--threads'] = int(arguments["--threads"])

	if arguments['--cluster']:
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

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
