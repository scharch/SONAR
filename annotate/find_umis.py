#!/usr/bin/env python3

"""
find_umis.py

This is a helper script to split up UMI identification from large sequencing runs
    for the sake of speed and memory usage.

Usage: find_umis.py FASTA FORMAT [ --cell 0,16 --umi 16,26 --r2umi 0,8 ] [ --pe --revcomp ] [ --cellWhiteList barcodes.txt | --cellPattern NNNNNN ] [ --umiWhiteList barcodes.txt | --umiPattern NNNNNN ] [ --umi2WhiteList barcodes.txt | --umi2Pattern NNNNNN ] [ --minQ Q ]

Options:
    FASTA                          Subsampled fasta/q file produced by 1.0-preprocess.py
    FORMAT                         Format of input file, fasta or fastq
    --cell 0,16                    See 1.0-preprocess.py for explanation
    --umi 16,26                    See 1.0-preprocess.py for explanation
    --r2umi 0,8                    See 1.0-preprocess.py for explanation
    --pe                           Flag to indicate use of 10x PE short reads strategy for feature barcoding.
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

		count     = 0
		bad_umi   = 0
		low_qual  = 0
		print("%s: Starting to look for UMIs in %s" % (datetime.datetime.now(), arguments["FASTA"]) )

		if arguments['--pe']:
			r2Handle = open( "r2" + arguments["FASTA"], "r" )
			r2Parser = SeqIO.parse( r2handle, arguments["FORMAT"])

		with open(arguments["FASTA"], "r") as handle:
			for seq in SeqIO.parse( handle, arguments["FORMAT"]):
				count += 1

				#check for derep-ed ness
				reads = 1
				#check_derep = re.search(";size=(\d+)", seq.id)
				#if check_derep:
				#	reads = int( check_derep.group(1) )
				#	seq.id = re.sub(";size=(\d+)", "", seq.id) #I don't remember why this is here, seems like not the correct behavior

				cell_barcode = str(seq.seq[ cb_start:cb_end ])
				fwd_id       = str(seq.seq[ umi_start:umi_end ])
				rev_id       = str(seq.seq.reverse_complement()[ umi2_start:umi2_end ])

				#check whitelists/patterns
				if cell_barcode != "":
					if arguments["FORMAT"]=="fastq" and any([ x<arguments['--minQ'] for x in seq.letter_annotations['phred_quality'][cb_start:cb_end] ]):
						low_qual += 1
						continue
					elif arguments['--cellWhiteList'] is not None:
						if not cell_barcode in cellWhiteList:
							bad_umi += 1
							continue
					elif arguments['--cellPattern'] is not None:
						if not re.match(arguments['--cellPattern'], cell_barcode):
							bad_umi += 1
							continue

				if fwd_id != "":
					if arguments["FORMAT"]=="fastq" and any([ x<arguments['--minQ'] for x in seq.letter_annotations['phred_quality'][umi_start:umi_end] ]):
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
					if arguments["FORMAT"]=="fastq" and any([ x<arguments['--minQ'] for x in seq.reverse_complement().letter_annotations['phred_quality'][umi2_start:umi2_end] ]):
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
				seq = seq[ max(cb_end, umi_end): ]
				if umi2_end > 0:
					seq = seq[ : -umi2_end]
				if arguments['--pe']:
					tempseq = next(r2Parser)
					if re.sub("/1$","",seq.id) == re.sub("/2$","", tempseq.id):
						seq = tempseq
					else:
						sys.exit( f"Error: sequence id mismatch between R1 and R2 in {arguments['FASTA']}: {seq.id} vs {tempseq.id}" )
				elif arguments['--revcomp']:
					seq.seq = seq.seq.reverse_complement()

				if cell_barcode != "":
					seq.id += ";cell=%s"%cell_barcode
				if molecule_id	!= "":
					seq.id += ";umi=%s"%molecule_id
				else:
					#no umi, but store cell barcode as the umi in the data structure to prevent errors
					molecule_id = cell_barcode

				if (cell_barcode, molecule_id) not in umi_dict:
					umi_dict[ (cell_barcode, molecule_id) ] = { 'cell':cell_barcode, 'umi':molecule_id, 'count':reads, 'seqs':[seq] }
				else:
					umi_dict[ (cell_barcode, molecule_id) ]['count'] += reads
					umi_dict[ (cell_barcode, molecule_id) ]['seqs'].append(seq)

			print( "%s: Finished %s: %d sequences in %d UMIs; Discarded %d reads with low quality UMIs and %d additional reads with illegal UMIs." % (datetime.datetime.now(), arguments["FASTA"], count, len(umi_dict), low_qual, bad_umi) )

			with open(re.sub(arguments["FORMAT"],"pickle",arguments["FASTA"]), 'wb') as pickle_out:
				pickle.dump( umi_dict, pickle_out )



if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--minQ']     = int( arguments['--minQ'] )

	prj_tree = ProjectFolders( os.getcwd() )
	prj_name = fullpath2last_folder(prj_tree.home)

	#log command line
	#logCmdLine(sys.argv)

	iupac = { "A":"A", "C":"C", "G":"G", "T":"[UT]", "U":"[UT]", "M":"[AC]", "R":"[AG]", "W":"[AT]", "S":"[CG]", "Y":"[CT]", "K":"[GT]", "V":"[ACG]", "H":"[ACT]", "D":"[AGT]", "B":"[CGT]", "N":"[ACGTU]" }
	cellWhiteList = []
	umiWhiteList  = []
	umi2WhiteList = []
	if arguments['--cellWhiteList'] is not None:
		with open(arguments['--cellWhiteList'], "r") as codes:
			for bc in codes.readlines():
				cellWhiteList.append(bc.strip())
	elif arguments['--cellPattern'] is not None:
		arguments['--cellPattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--cellPattern'])

	if arguments['--umiWhiteList'] is not None:
		with open(arguments['--umiWhiteList'], "r") as codes:
			for bc in codes.readlines():
				umiWhiteList.append(bc.strip())
	elif arguments['--umiPattern'] is not None:
		arguments['--umiPattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--umiPattern'])

	if arguments['--umi2WhiteList'] is not None:
		with open(arguments['--umi2WhiteList'], "r") as codes:
			for bc in codes.readlines():
				umi2WhiteList.append(bc.strip())
	elif arguments['--umi2Pattern'] is not None:
		arguments['--umi2Pattern'] = re.sub("\w", lambda x: iupac[x.group().upper()], arguments['--umi2Pattern'])


	main()
