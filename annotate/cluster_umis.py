#!/usr/bin/env python3

"""
cluster_umis.py

This is a helper script to split up UMI consensus generation.

Usage: cluster_umis.py PICKLE MINSIZE DIR [ --isCell --isFeature ]

Options:
    PICKLE          Pickled umi dictionary produced by 1.0-preprocess.py
    MINSIZE         Minimum reads/umi or umis/cell
    DIR             Directory in which to carry out clustering
    --isCell        Flag to indicate generation of cell metaconsenus (instead of
	                    individual UMI consensus) [default: False]
    --isFeature     Flag to adjust clustering for short feature barcoding oligos
	                    instead of V(D)Js [default: False]

Split out from 1.0-preprocess.py by Chaim A Schramm on 2019-06-18.
Added option for feature barcodes by CA Schramm 2019-10-08.
Tweaked clustering thresholds by CAS 2020-02-04.
Changed minUMIs to a per metaconsenus (instead of per cell) threshold by CAS 2020-02-12.
Changed structure of umi_dict by CAS 2020-08-07.
Added a level to the subdirectory structure and fixed a minor bug that was
    reporting phantom cells having "fewer than 1 UMI" by CAS 2020-12-17.


Copyright (c) 2019-2020 Vaccine Research Center, National Institutes of Health, USA.
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
 		umi_iter = pickle.load(pickle_in)

	results = {}
	small	= 0
	multi	= 0

	print( f"{datetime.datetime.now()}: Generating consensus sequences from {arguments['PICKLE']}..." )

	for cell, cell_dict in umi_iter:

		#kludge to make data structure the same
		if arguments['--isCell']:
			cell_dict = { cell:cell_dict }

		for umi in cell_dict:

			#for isCell==True, numReads is actually 'numUMIs'
			numReads = sum( [ len(cell_dict[umi][s]) for s in cell_dict[umi] ] )
			numSeqs  = len( cell_dict[umi]	)

			#check thresholds
			if numReads < arguments['MINSIZE']:
				small += 1  #not entirely correct for metaconsenus case because it might be 
							#     more than one metaconsensus that doesn't have enough UMIs
				continue

			if numSeqs == 1:
				#save time on singletons (if they weren't excluded by the read threshold)
				seq   = next(iter(cell_dict[umi]))
				seqid = cell_dict[umi][seq][0]

				if arguments['--isCell']:
					totalReads = sum( [ int(cc.group(1)) for cc in [ re.search("consensus_count=(\d+)",mi) for mi in cell_dict[umi][seq] ] ] )
					seqid = "%s.1 cell_id=%s duplicate_count=%d consensus_count=%d"%( cell, cell, numReads, totalReads )
				else:
					seqid += ";seqs=1;size=%d;consensus_count=%d" % (numReads, numReads)

				if not cell in results:
					results[cell] = defaultdict( list )
				results[ cell ][ seq ].append( seqid )

			else:

				#use the first half of the cell barcode as an intermediate
				#directory to avoid breaking the file system with large
				#datasets - the number of "cells" that we may have to process
				#here is many fold greater than the 10s of 1000s of actual
				#cells in an experiment
				intermed = cell[ : int(len(cell)/2) ]
				subdir = arguments['DIR'] + "/" + intermed

				#iddef  = "3"
				if not arguments['--isCell']:
					subdir += "/%s" % cell
					#iddef = "2"

				#cluster and rapid align with vsearch
				os.makedirs(subdir, exist_ok=True)
				with open("%s/%s.fa" % (subdir, umi), "w") as handle:
					for s in cell_dict[umi]:
						if arguments['--isCell']:
							#second round clustering, don't collapse identical sequences
							for mi in cell_dict[umi][s]:
								handle.write(f">{mi}\n{s}\n")
						else:
							#first round clustering, only write once per sequence
							#   but need to add size annotation
							handle.write(f">{cell_dict[umi][s][0]};size={len(cell_dict[umi][s])}\n{s}\n")

				otherOpts = [ "-mincols", '150' ] #VDJ alignment overlap
				if arguments['--isFeature']:
					otherOpts = [ "-minseqlength", '10' ] #in case of naked feature barcodes
				if arguments['--isCell']:
					otherOpts += ['-id', '0.99', '-uc', f"{subdir}/{umi}.uc"]
				else:
					otherOpts += ['-id', '0.95', '-sizein'] #don't want to weight by abundance for cell metaconsensus

				subprocess.call([vsearch,
						 "-cluster_fast", "%s/%s.fa" % (subdir, umi),
						 "-consout", "%s/%s_cons.fa" % (subdir, umi),
						 "-iddef", "3",
						 "-sizeout",
						 "-gapopen", "10I/10E", #lower gap open penalty to better account for internal indels
						 "-gapext", "2I/2E", #don't make endgaps cheaper; encourages TSOs to align properly
						 "-clusterout_sort", #so we can look at just the biggest
						 "-quiet" #supress screen clutter
						 ] + otherOpts )

				if arguments['--isCell']:
					totalReads = dict()
					with open(f"{subdir}/{umi}.uc", 'r') as ucfile:
						reader = csv.reader(ucfile, delimiter="\t")
						for row in reader:
							if row[0]=="C": break
							tr = re.search("size=(\d+)", row[8])
							if tr:
								if row[9] == "*":
									totalReads[row[8].split(";")[0]] = int(tr.group(1))
								else:
									totalReads[row[9].split(";")[0]] += int(tr.group(1))

				with open("%s/%s_cons.fa" % (subdir, umi), 'r') as cons_file:
					seq_number = 0
					for cons in SeqIO.parse(cons_file, "fasta"):
						seq_number += 1
						cons.id = re.sub("centroid=","",cons.id,count=1)
						if not arguments['--isCell'] and seq_number > 1:
							#how to handle more than one cluster per umi?
							#  -depends on presence/absence of cell barcodes, I guess. user param?
							#use it to do error checking???
							multi += 1
							break

						clusterReads = re.search(";seqs=(\d+);size=(\d+)",cons.id)
						if clusterReads:
							if int(clusterReads.group(2)) < arguments['MINSIZE']:
								small += 1
								continue

							if arguments['--isCell']:
								cons.id	 = "%s.%d cell_id=%s duplicate_count=%s consensus_count=%s"%( cell, seq_number, cell, clusterReads.group(1), totalReads[cons.id.split(";")[0]] )
							else:
								cons.id += ";consensus_count=%s" % clusterReads.group(2) #save size annotation for further clustering/dereplication

							cons.description = ""
							seq = str(cons.seq)
							if not cell in results:
								results[cell] = defaultdict( list )
							results[ cell ][ seq ].append( cons.id )

	with open(re.sub("cons_in","cons_out",arguments["PICKLE"]), 'wb') as pickle_out:
		pickle.dump( {'results':results, 'small':small, 'multi':multi}, pickle_out )

	print( f"{datetime.datetime.now()}: Finished processing {arguments['PICKLE']}" )


if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['MINSIZE'] = int( arguments['MINSIZE'] )

	prj_tree = ProjectFolders( os.getcwd() )
	prj_name = fullpath2last_folder(prj_tree.home)

	#log command line
	#logCmdLine(sys.argv)


	main()
