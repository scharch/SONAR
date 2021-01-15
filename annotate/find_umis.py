#!/usr/bin/env python3

"""
find_umis.py

This is a helper script to split up UMI identification from large sequencing runs
    for the sake of speed and memory usage.

Usage: find_umis.py FASTA CHUNKSIZE CHUNKNUM OUTPUT [ --cell 0,16 --umi 16,26 --r2umi 0,8 ] [ --pe r2.fq | --revcomp ] [ --cellWhiteList barcodes.txt | --cellPattern NNNNNN ] [ --umiWhiteList barcodes.txt | --umiPattern NNNNNN ] [ --umi2WhiteList barcodes.txt | --umi2Pattern NNNNNN ] [ --minQ Q ]

Options:
    FASTA                          Input fasta/q file
    CHUNKSIZE                      How many sequence to process in this worker.
    CHUNKNUM                       Which number worker is this
    OUTPUT                         File name for saving the pickled output.
    --cell 0,16                    See 1.0-preprocess.py for explanation
    --umi 16,26                    See 1.0-preprocess.py for explanation
    --r2umi 0,8                    See 1.0-preprocess.py for explanation
    --pe r2.fq                     File name of R2 for 10x PE short reads strategy for feature barcoding.
    --revcomp                      Flag to reverse-complement the sequence after UMI identification for SE feature barcoding.
    --cellWhiteList barcodes.txt   See 1.0-preprocess.py for explanation
    --cellPattern NNNNNN           See 1.0-preprocess.py for explanation
    --umiWhiteList barcodes.txt    See 1.0-preprocess.py for explanation
    --umiPattern NNNNNN            See 1.0-preprocess.py for explanation
    --umi2WhiteList barcodes.txt   See 1.0-preprocess.py for explanation
    --umi2Pattern NNNNNN           See 1.0-preprocess.py for explanation
    --minQ Q                       See 1.0-preprocess.py for explanation

Split out from 1.0-preprocess.py by Chaim A Schramm on 2019-06-18.
Added PE and REVCOMP flags for handling feature barcoding by CA Schramm 2019-10-08.
Changed structure of umi_dict by CAS 2020-08-07.
Changed whitelists sets instead of lists by CAS 2020-08-07.
Using pyfastx indexing instead of pre-splitting input files by CAS 2020-12-16.

Copyright (c) 2019-2020 Vaccine Research Center, National Institutes of Health, USA.
    All rights reserved.

"""

import sys, os, re, pickle, gzip
from functools import partial
from docopt import docopt
import datetime
from collections import defaultdict
from Bio import SeqIO
import pyfastx
from numpy import random
from time import sleep

try:
    from SONAR.annotate import *
except ImportError:
    find_SONAR = sys.argv[0].split("SONAR/annotate")
    sys.path.append(find_SONAR[0])
    from SONAR.annotate import *

def main():

		cb_start, cb_end = 0, 0
		if arguments['--cell'] is not None:
			cb_start, cb_end = [ int(x) for x in arguments['--cell'].split(",") ]

		umi_start, umi_end = 0, 0
		if arguments['--umi'] is not None:
			umi_start, umi_end = [ int(x) for x in arguments['--umi'].split(",") ]

		umi2_start, umi2_end = 0, 0
		if arguments['--r2umi'] is not None:
			umi2_start, umi2_end = [ int(x) for x in arguments['--umi2'].split(",") ]

		#start reading in the files
		umi_dict    = {}
		bad_umi   = 0
		low_qual  = 0

		# #introduce a small random delay to try and prevent segfaults
		# #when loading the index file
		# sleeptime = 15 * random.uniform(0, 10)
		# sleep(sleeptime)

		try:
			seqFile = pyfastx.Fastq( arguments['FASTA'])
		except RuntimeError:
			seqFile = pyfastx.Fasta(arguments['FASTA'])

		if arguments['--pe'] is not None:
			try:
				r2File = pyfastx.Fastq( arguments['--pe'])
			except RuntimeError:
				r2File = pyfastx.Fasta(arguments['--pe'])

		print("loaded R2", flush=True)

		start = arguments['CHUNKNUM']*arguments['CHUNKSIZE']
		stop  = min( len(seqFile), start+arguments['CHUNKSIZE'] ) - 1

		for seqNum in range(start, stop):

			seq = seqFile[seqNum]

			cell_barcode = seq.seq[ cb_start:cb_end ]
			fwd_id       = seq.seq[ umi_start:umi_end ]
			#try:
			rev_id       = seq.antisense[ umi2_start:umi2_end ]
			# except:
			# 	#temporary kludge until pyfastx supports `antisense` on fastqs
			# 	rev_id = str( Seq.Seq(seq.seq).reverse_complement()[umi2_start:umi2_end] )

			#check whitelists/patterns
			if cell_barcode != "":
				if hasattr(seq, 'quali') and any([ x<arguments['--minQ'] for x in seq.quali[cb_start:cb_end] ]):
					low_qual += 1
					continue
				if arguments['--cellWhiteList'] is not None:
					if not cell_barcode in cellWhiteList:
						bad_umi += 1
						continue
				elif arguments['--cellPattern'] is not None:
					if not re.match(arguments['--cellPattern'], cell_barcode):
						bad_umi += 1
						continue

			if fwd_id != "":
				if hasattr(seq, 'quali') and any([ x<arguments['--minQ'] for x in seq.quali[umi_start:umi_end] ]):
					low_qual += 1
					continue
				elif arguments['--umiWhiteList'] is not None:
					if not fwd_id in umiWhiteList:
						bad_umi += 1
						continue
				elif arguments['--umiPattern'] is not None:
					if not re.match(arguments['--umiPattern'], fwd_id):
						bad_umi += 1
						continue

			if rev_id != "":
				if hasattr(seq, 'quali') and any([ x<arguments['--minQ'] for x in seq.quali[-umi2_end:-umi2_start] ]):
					low_qual += 1
					continue
				elif arguments['--umi2WhiteList'] is not None:
					if not rev_id in umi2WhiteList:
						bad_umi += 1
						continue
				elif arguments['umi2--Pattern'] is not None:
					if not re.match(arguments['--umi2Pattern'], rev_id):
						bad_umi += 1
						continue

			#combine UMIs and trim them from sequence
			molecule_id = fwd_id + rev_id

			if arguments['--pe']:
				finalSeq = r2File[seqNum].seq
			else:
				leftEnd  = max(cb_end, umi_end)
				rightEnd = None
				if umi2_end > 0:
					rightEnd = umi2_end
				if arguments['--revcomp']:
					if rightEnd == None:
						rightEnd = 0
					#try:
					finalSeq = seq.antisense[ -rightEnd : -leftEnd ]
					# except:
					# 	#temporary kludge until pyfastx supports `antisense` on fastqs
					# 	finalSeq = str( Seq.Seq(seq.seq).reverse_complement()[ -rightEnd : -leftEnd ] )
				else:
					finalSeq = seq.seq[ leftEnd : rightEnd ]

			finalID = seq.name
			if cell_barcode != "":
				finalID += ";cell=%s"%cell_barcode
			else:
				#no cell barcode, but store umi as the cell barcode in the data structure to prevent errors
				cell_barcode = molecule_id
				
			if molecule_id	!= "":
				finalID += ";umi=%s"%molecule_id
			else:
				#no umi, but store cell barcode as the umi in the data structure to prevent errors
				molecule_id = cell_barcode

			if cell_barcode not in umi_dict:
				umi_dict[cell_barcode] = defaultdict( dict )

			if finalSeq not in umi_dict[cell_barcode][molecule_id]:
				umi_dict[ cell_barcode ][ molecule_id ][ finalSeq ] = [ finalID ]
			else:
				umi_dict[ cell_barcode ][ molecule_id ][ finalSeq ].append(finalID)

		print( "%s: Finished %s: %d sequences in %d UMIs; Discarded %d reads with low quality UMIs and %d additional reads with illegal UMIs." % (datetime.datetime.now(), arguments["FASTA"], arguments['CHUNKSIZE'], len(umi_dict), low_qual, bad_umi) )

		with open( arguments["OUTPUT"], 'wb') as pickle_out:
			pickle.dump( umi_dict, pickle_out )



if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['CHUNKSIZE'] = int( arguments['CHUNKSIZE'] )
	arguments['CHUNKNUM']  = int( arguments['CHUNKNUM'] )
	arguments['--minQ']    = int( arguments['--minQ'] )

	prj_tree = ProjectFolders( os.getcwd() )
	prj_name = fullpath2last_folder(prj_tree.home)

	#log command line
	#logCmdLine(sys.argv)

	iupac = { "A":"A", "C":"C", "G":"G", "T":"[UT]", "U":"[UT]", "M":"[AC]", "R":"[AG]", "W":"[AT]", "S":"[CG]", "Y":"[CT]", "K":"[GT]", "V":"[ACG]", "H":"[ACT]", "D":"[AGT]", "B":"[CGT]", "N":"[ACGTU]" }
	cellWhiteList = set()
	umiWhiteList  = set()
	umi2WhiteList = set()
	if arguments['--cellWhiteList'] is not None:
		with open(arguments['--cellWhiteList'], "r") as codes:
			for bc in codes.readlines():
				cellWhiteList.add(bc.strip())
	elif arguments['--cellPattern'] is not None:
		arguments['--cellPattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--cellPattern'])

	if arguments['--umiWhiteList'] is not None:
		with open(arguments['--umiWhiteList'], "r") as codes:
			for bc in codes.readlines():
				umiWhiteList.add(bc.strip())
	elif arguments['--umiPattern'] is not None:
		arguments['--umiPattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--umiPattern'])

	if arguments['--umi2WhiteList'] is not None:
		with open(arguments['--umi2WhiteList'], "r") as codes:
			for bc in codes.readlines():
				umi2WhiteList.add(bc.strip())
	elif arguments['--umi2Pattern'] is not None:
		arguments['--umi2Pattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--umi2Pattern'])


	main()
