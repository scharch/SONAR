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

Copyright (c) 2019 Vaccine Research Center, National Institutes of Health, USA.
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

	for umi in umi_iter:

		#check read threshold
		if umi['count'] < arguments['MINSIZE']:
			small += 1
			continue

		if len(umi['seqs']) == 1:
			#save time on singletons (if they weren't excluded by the read threshold)
			if arguments['--isCell']:
				umi['seqs'][0].id = "%s.1 cell_id=%s duplicate_count=1 consensus_count=%s"%( umi['cell'], umi['cell'], umi['count'] )
				umi['seqs'][0].description = ""
			else:
				umi['seqs'][0].id += ";seqs=1;size=%d;consensus_count=%d" % (umi['count'],umi['count'])
				umi['seqs'][0].description = ""

			if (umi['cell']) not in results:
				results[ umi['cell'] ] = { 'cell':umi['cell'], 'umi':umi['cell'], 'count':1, 'seqs':umi['seqs'].copy() }
			else:
				results[ umi['cell'] ]['count'] += 1
				results[ umi['cell'] ]['seqs']  += umi['seqs']

		else:

			subdir = arguments['DIR']
			#iddef  = "3"
			if not arguments['--isCell']:
				subdir += "/%s" % umi['cell']
				#iddef = "2"

			#cluster and rapid align with vsearch
			os.makedirs(subdir, exist_ok=True)
			with open("%s/%s.fa" % (subdir, umi['umi']), "w") as handle:
				SeqIO.write(umi['seqs'], handle, "fasta")

			lenOpts = [ "-mincols", '150' ] #VDJ alignment overlap
			if arguments['--isFeature']:
				lenOpts = [ "-minseqlength", '10' ] #in case of naked feature barcodes

			subprocess.call([vsearch,
					 "-cluster_fast", "%s/%s.fa" % (subdir, umi['umi']),
					 "-consout", "%s/%s_cons.fa" % (subdir, umi['umi']),
					 "-id", "0.97",
					 "-iddef", "3",
					 "-sizein", "-sizeout",
					 "-gapopen", "10I/10E", #lower gap open penalty to better account for internal indels
					 "-gapext", "2I/2E", #don't make endgaps cheaper; encourages TSOs to align properly
					 "-clusterout_sort", #so we can look at just the biggest
					 "-quiet" #supress screen clutter
					 ] + lenOpts )

			with open("%s/%s_cons.fa" % (subdir, umi['umi']), 'r') as cons_file:
				seq_number = 0
				for cons in SeqIO.parse(cons_file, "fasta"):
					seq_number += 1
					if not arguments['--isCell'] and seq_number > 1:
						#how to handle more than one cluster per umi?
						#  -depends on presence/absence of cell barcodes, I guess. user param?
						#use it to do error checking???
						multi += 1
						break

					num_reads = re.search(";seqs=(\d+);size=(\d+)",cons.id)
					if num_reads:
						if arguments['--isCell']:
							cons.id	 = "%s.%d cell_id=%s duplicate_count=%s consensus_count=%s"%( umi['cell'], seq_number, umi['cell'], num_reads.group(1), num_reads.group(2) )
						else:
							if int(num_reads.group(2)) < arguments['MINSIZE']:
								small += 1
								continue
							else:
								cons.id += ";consensus_count=%s" % num_reads.group(2) #save size annotation for further clustering/dereplication

						cons.description = ""
						if (umi['cell']) not in results:
							results[ umi['cell'] ] = { 'cell':umi['cell'], 'umi':umi['cell'], 'count':1, 'seqs':[cons] }
						else:
							results[ umi['cell'] ]['count'] += 1
							results[ umi['cell'] ]['seqs'].append(cons)

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
