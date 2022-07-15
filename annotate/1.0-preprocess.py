#!/usr/bin/env python3

"""
1.0-preprocess.py

This script preprocesses raw sequencing data in preparation for VDJ analysis. Paired end fastq
    files will be filtered and merged using vsearch, or single read fastq or fasta files can be
    specified. For fastq input, an optional QC analysis can be printed (from vsearch --fastq_stats)
    by using the --printQC option.
After QC, this script proceeds to look for cell barcodes and UMIs. Currently, cell barcodes are
    only supported on the 5' end of the read; UMIs can be on either end, including support for
    double UMI protocols. (Note that it treats both ends equivalently, so it's not quite compatible
    with, eg, the Reddy lab's MAF protocol.) The script allows for filtering barcodes/UMIs against
    either a known whitelist or a designed pattern, but does not currently check for UMI
    error/splitting. (I plan to add it, but it's mostly moot in the context of a cell barcode.)
The script then uses vsearch to generate a consensus sequence for each UMI and, importantly, a set
    of meta-consensus sequences for the UMIs in each cell. In order to be conservative, this will
    sometimes (frequently) generate multiple copies of the same Ig transcript in the same cell. Use
    1.5-single_cell_statistics.py after annotation to clean things up a bit.

TODO:
   * Add dereplication and UMI collision detection when no cell barcodes are present
   * Allow UMI reuse under certain conditions?
   * Should there be a threshold (10:1? 100:1?) at which I ignore conflicting sample assignments
         from cell hashing and just go with the majority?
   * Come up with a better way to decide which cells are "real"; 1 UMI might be okay if it has 50
         reads; OTOH 3 UMIs might not be enough if they do not all cluster together.

Usage: 1.0-preprocess.py --input read1.fq... [ --reverse read2.fq... ] [ ( --featureLibrary features.fq... --featureList features.tsv ) ] [ --cellWhiteList barcodes.txt | --cellPattern NNNNNN ] [ --umiWhiteList barcodes.txt | --umiPattern NNNNNN ] [ --umi2WhiteList barcodes.txt | --umi2Pattern NNNNNN ] [ options ]

Options:
    --input read1.fq               File with raw data to process. Can be used multiple times **IF**
                                       files represent technical or sequencing replicates, where
                                       overlapping cell barcodes/UMIs can be assumed to be derived
                                       from the same physical source.
    --reverse read2.fq             File with reverse reads for paired ends. If used, must be provided
                                       the same number of times (and in the same order) as --input.
    --featureLibrary features.fq   File with sequence reads from cell hashing/feature barcoding
                                       libraries. Can be specified twice for 10x-recommended 26+25
                                       seqeuncing strategies, in which case the cell barcode and UMI
                                       will be looked for on R1 and the feature barcodes on R2. If
                                       only one file is input, the script will look for the reverse
                                       complements of the feature barcodes downstream of the UMI.
    --cell 0,16                    Python-style zero-indexed, semi-open interval with the expected
                                       position of the cell barcode, if present.
    --umi 16,26                    Python-style zero-indexed, semi-open interval with the expected
                                       position of the UMI, if present.
    --r2umi 0,8                    Python-style zero-indexed, semi-open interval with the expected
                                       position of the UMI on Read 2, if present. If only one file
                                       is input, the script will treat the 3' end of the read as a
                                       reverse-complemented R2 and look for the umi there.
    --cellWhiteList barcodes.txt   A file with allowed cell barcodes, one per line. Reads with
                                       apparent barcodes that are not on the list will be discarded.
                                       Mutually exclusive with cellPattern.
    --cellPattern NNNNNN           A nucleotide pattern (using IUPAC ambiguity codes) describing the
                                       cell barcodes. Reads with apparent barcodes that do not match
                                       the pattern will be discarded. Mutually exclusive with
                                       cellWhiteList.
    --umiWhiteList barcodes.txt    A file with allowed UMIs, one per line. Reads with apparent UMIs
                                       that are not on the list will be discarded. Mutually
                                       exclusive with umiPattern.
    --umiPattern NNNNNN            A nucleotide pattern (using IUPAC ambiguity codes) describing the
                                       UMIs. Reads with apparent UMIs that do not match the pattern
                                       will be discarded. Mutually exclusive with umiWhiteList.
    --umi2WhiteList barcodes.txt   A file with allowed R2 UMIs, one per line. Reads with apparent R2
                                       UMIs that are not on the list will be discarded. Mutually
                                       exclusive with umi2Pattern.
    --umi2Pattern NNNNNN           A nucleotide pattern (using IUPAC ambiguity codes) describing the
                                       R2 UMIs. Reads with apparent R2 UMIs that do not match the
                                       pattern will be discarded. Mutually exclusive with
                                       umi2WhiteList.
    --featureList features.tsv     A tab-delimited text file with cell hashing/feature-barcoding
                                       oligos in the first column and their respective features in
                                       the second column. For cell hashing specifically, use
                                       "sample:pool1," which will produce a "sample" column in the
                                       rearrangments.tsv; all other values will be used as the name of
                                       a custom column in the cell_stats.tsv with the value
                                       corresponding to the number of detected UMIs.
    --filterOptions options        A string of options to be passed to vsearch --fastx_filter for
                                       quality control Will be applied equally to R1 and R2 *before*
                                       merging, or to a single input file if no R2 is specified.
                                       [default: -fastq_truncee 5]
    --mergeOptions options         A string of options to be passed to vsearch --fastq_mergepairs.
                                       Ignored if only one input file is specified.
                                       [default: -fastq_minmergelen 350 -fastq_maxdiffs 100 -fastq_maxdiffpct 25 -fastq_eeout -fasta_width 0]
    --printQC file.log             A file in which to save a report on the quality of the final
                                       input sequences, after QC but before UMI processing,
                                       using vsearch --fastq_stats.
    --logFile preprocess.log       Where to save the log file. [default: output/logs/preprocess.log]
    --minQ 20                      Minimum PHRED score for all bases in a UMI or cell barcode. Reads
                                       with *any* base in the UMI/barcode below this threshold will
                                       be discarded. [default: 20]
    --minReads 1                   Minimum number of reads per UMI. UMIs with fewer reads will be
                                       discarded. [default: 3]
    --minUMIs 1                    Minimum number of UMIs per metaconsenus/final sequence. In theory,
                                       a value >1 should help remove background contamination. [default: 1]
    --umiOutput file.fa            File in which to save the UMI-processed sequences. Ignored if
                                       UMIs are not present. [default: byUMI.fa]
    --cellOutput file.fa           File in which to save the cell-processed sequences. Ignored if
                                       cell barcodes are not present. [default: byCell.fa]
    --runVBlast                    Flag to call 1.1 script when finished. Additional options to that
                                       script (including the ability to automatically call 1.2 and
                                       further scripts in Module 1) are listed below. This script
                                       will not check the validity of options passed downstream, so
                                       user beware.
    --cluster                      Flag to submit chunk jobs to cluster instead of running them
                                       locally. [default: False]
    --threads 1                    Number of threads to use. [default: 1]
    -f                             Flag to force overwriting of old files. [default: False]
    --keepWorkFiles                Flag to prevent deletion of intermediate files (useful for
                                       inspecting UMI clustering). [default: False]

Options for other annotation scripts (see those help messages for details):
    --locus H
    --species human
    --lib LIB
    --derep
    --minl <300>
    --maxl <600>
    --npf <10000>
    --runJBlast
    --jlib LIB
    --dlib LIB
    --clib LIB
    --noD
    --noC
    --runFinalize
    --jmotif TGGGG
    --nterm OPT
    --noclean
    --noFallBack
    --runClustering
    --file FILE
    --min1 1
    --min2 3
    --id .99
    --maxgaps 0
    --runCellStatistics
    --rearrangements rearrangements.tsv
    --save OPT

Created by Chaim A Schramm on 2019-01-09.
Code refactored and many options added by CAS 2019-02-12.
Added options for fasta input and/or multiple input files by CAS 2019-03-04.
Updated how Module 1 scripts chain together by CA Schramm 2019-04-01.
Fixed bug for singletons by CA Schramm 2019-05-02.
Switched to `cluster_fast` and changed gap penalties by CAS 2019-05-09.
Changed clustering threshold to 97% by CA Schramm 2019-05-22.
Fixed bug that was ignoring `--minQ` by CAS 2019-06-18.
Added `--keepWorkFiles` flag by CAS 2019-06-18.
Split off `find_umis.py` and `cluster_umis.py` into separate helper scripts
    and option to parallelize these on a cluster by CAS 2019-06-19.
Added code for feature barcoding by CAS 2019-10-08.
Added species option by CAS 2020-02-06.
Changed minUMIs to a per metaconsenus threshold and set default to 1 (for now)
    by CAS 2020-02-12.

Copyright (c) 2019-2020 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, shutil, gzip, csv
from docopt import docopt
import itertools
import pickle
from collections import defaultdict
from Bio import SeqIO
from multiprocessing import Pool
import datetime
from functools import partial
import random, math
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


def callFinder(num, format, args, stem="chunk"):
	cmd = f"{SCRIPT_FOLDER}/annotate/find_umis.py  {prj_tree.preprocess}/{stem}{num:04}.{format} {format} {args}"
	os.system( cmd )


def getUmiConsensus(num, minSize, workdir, isCell=False, feature=False, clustType="umi"):
	if isCell:
		clustType = "cell"
	cmd = f"{SCRIPT_FOLDER}/annotate/cluster_umis.py {prj_tree.preprocess}/{clustType}_cons_in_{num:04}.pickle {minSize} {workdir}"
	if isCell:
		cmd += " --isCell"
	elif feature:
		cmd += " --isFeature"
	os.system( cmd )


def processFeatures():

	fileformat = "fastq"
	try:
		#gzip?
		if re.search("gz$", arguments['--featureLibrary'][0]):
			_open = partial(gzip.open,mode='rt')
		else:
			_open = partial(open, mode='r')
		with _open(arguments['--featureLibrary'][0]) as checkInput:
			parser = SeqIO.parse(checkInput, "fastq")
			testSeq = next(parser)
	except StopIteration:
		#fasta input
		fileformat		     = "fasta"


	fInd = 0
	with _open(arguments['--featureLibrary'][0]) as toSplit:
		splitForID = SeqIO.parse( toSplit, format=fileformat )
		for chunk in iterator_slice(splitForID, 50000):
			fInd += 1
			with open( "%s/features%04d.%s"%(prj_tree.preprocess, fInd, fileformat), 'w' ) as writeChunk:
				SeqIO.write( chunk, writeChunk, fileformat )
	#10x PE short read strategy
	if len(arguments['--featureLibrary']) == 2:
		ind2 = 0
		with _open(arguments['--featureLibrary'][1]) as toSplit:
			splitForID = SeqIO.parse( toSplit, format=fileformat )
			for chunk in iterator_slice(splitForID, 50000):
				ind2 += 1
				with open( "%s/r2features%04d.%s"%(prj_tree.preprocess, fInd, fileformat), 'w' ) as writeChunk:
					SeqIO.write( chunk, writeChunk, fileformat )


	#construct the command for find_umis
	featureOpts = ""
	for opt in ['--cell', '--umi', '--r2umi', '--cellWhiteList', '--cellPattern', '--umiWhiteList', '--umiPattern', '--umi2WhiteList', '--umi2Pattern', '--minQ' ]:
		if arguments[opt] is not None:
			featureOpts += " %s '%s'" % (opt, arguments[opt])
	if len(arguments['--featureLibrary']) == 2:
		featureOpts += " --pe"
	else:
		featureOpts += " --revcomp"

	#call find_umis either on cluster or locally
	if arguments['--cluster']:
		with open("%s/featuresjob.sh"%prj_tree.preprocess, 'w') as jobHandle:
			jobHandle.write(f"#!/bin/bash\n#$ -N featureUMIs\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/find_umis.py {prj_tree.preprocess}/features$NUM.{fileformat} {fileformat} {featureOpts}\n\n")
		subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%fInd, "%s/featuresjob.sh"%prj_tree.preprocess])
	else:
		partial_finder = partial( callFinder, format=fileformat, args=featureOpts, stem="features" )
		pool = Pool(arguments['--threads'])
		pool.map( partial_finder, range(1,fInd+1) )
		pool.close()
		pool.join()

	#collect output of find_umis
	feature_dict = {}
	pickleFiles = glob.glob(f"{prj_tree.preprocess}/feature*.pickle")
	for p in pickleFiles:
		with open(p, 'rb') as pickle_in:
			chunk_dict = pickle.load(pickle_in)
			for (cb, mi) in chunk_dict:
				if (cb, mi) not in feature_dict:
					feature_dict[ (cb, mi) ] = chunk_dict[ (cb, mi) ].copy()
				else:
					feature_dict[ (cb, mi) ]['count'] += chunk_dict[ (cb, mi) ]['count']
					feature_dict[ (cb, mi) ]['seqs']  += chunk_dict[ (cb, mi) ]['seqs']

	#generate pickles to pass to consensus algorithm
	dInd = 0
	for chunk in iterator_slice(feature_dict.values(), 5000):
		dInd += 1
		with open( f"{prj_tree.preprocess}/feature_cons_in_{dInd:04}.pickle", 'wb') as pickle_out:
			pickle.dump( chunk, pickle_out )

	#delete feature_dict to save memory
	feature_dict = None

	#spawn subprocesses
	if arguments['--cluster']:
		with open("%s/featurecons.sh"%prj_tree.preprocess, 'w') as jobHandle:
			jobHandle.write(f"#!/bin/bash\n#$ -N clusterFeatureUMIs\n#$-l h_vmem=32G\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/cluster_umis.py {prj_tree.preprocess}/feature_cons_in_$NUM.pickle {arguments['--minReads']} {prj_tree.preprocess} --isFeature\n\n")
		subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%dInd, "%s/featurecons.sh"%prj_tree.preprocess])
	else:
		partial_cons = partial( getUmiConsensus, minSize=arguments['--minReads'], workdir=prj_tree.preprocess, clustType="feature")

		pool = Pool(arguments['--threads'])
		blob = pool.map( partial_cons, range(1,dInd+1) )
		pool.close()
		pool.join()

	#read in feature barcode table
	hashingSeqs = dict()
	featureSeqs = dict()
	with open( arguments['--featureList'], 'r' ) as fHandle:
		reader = csv.reader(fHandle, delimiter="\t")
		for row in reader:
			if "N" in row[0]:
				sys.exit("I don't know how to handle a degenerate barcode, sorry!")
			checkHash = re.match("sample:(.*)", row[1])
			if checkHash:
				hashingSeqs[ row[0] ] = checkHash.group(1)
			else:
				featureSeqs[ row[0] ] = row[1]

	#collect output
	cellHashes   = defaultdict( dict )
	cellFeatures = defaultdict( dict )
	for p in glob.glob(f"{prj_tree.preprocess}/feature_cons_out_*.pickle"):
		with open(p, 'rb') as pickle_in:
			chunk_dict = pickle.load(pickle_in)
			for c in chunk_dict['results']:
				for s in chunk_dict['results'][c]['seqs']:
					for hash in hashingSeqs:
						if hash in str(s.seq):
							if hashingSeqs[hash] in cellHashes[c]:
								cellHashes[ c ][ hashingSeqs[hash] ] += 1
							else:
								cellHashes[ c ][ hashingSeqs[hash] ] = 1
							break

					for oligo in featureSeqs:
						if oligo in str(s.seq):
							if featureSeqs[oligo] in cellFeatures[c]:
								cellFeatures[ c ][ featureSeqs[oligo] ] += 1
							else:
								cellFeatures[ c ][ featureSeqs[oligo] ] = 1
							break

	if len(hashingSeqs) > 0:
		with open( f"{prj_tree.tables}/{prj_name}_hashes.tsv", 'w' ) as handle:
			writer = csv.writer( handle, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE )
			for cell in sorted(cellHashes.keys()):
				sample = "unknown"
				if len(cellHashes[cell]) > 1:
					sample = "ambiguous"
				elif len(cellHashes[cell]) == 1:
					sample = list(cellHashes[cell].keys())[0]
				writer.writerow( [cell, sample] )

	if len(featureSeqs) > 0:
		with open( f"{prj_tree.tables}/{prj_name}_features.tsv", 'w' ) as handle:
			writer = csv.writer( handle, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE)
			writer.writerow( ["cell_id"] + sorted(featureSeqs.values()) )
			for cell in sorted(cellFeatures.keys()):
				writer.writerow( [cell] + [ cellFeatures[cell].get(f, 0) for f in sorted(featureSeqs.values()) ] )

	print(str(datetime.datetime.now()) + " - Finished processing feature barcodes!")



def main():

	if len(arguments['--featureLibrary']) > 0:
		processFeatures()

	processedFiles = []

	for fileNum, inFile in enumerate(arguments['--input']):

		fileformat = "fastq"
		try:
			#gzip?
			if re.search("gz$", inFile):
				_open = partial(gzip.open,mode='rt')
			else:
				_open = partial(open, mode='r')
			with _open(inFile) as checkInput:
				parser = SeqIO.parse(checkInput, "fastq")
				testSeq = next(parser)
		except StopIteration:
			#fasta input (or maybe an empty file...)
			if len(arguments['--reverse']) > 0:
				sys.exit("Cannot merge PE reads in fasta format!")
			else:
				#ok, just don't try to do qual analysis
				arguments['--filterOptions'] = "None"
				arguments['--printQC']       =  None
				fileformat                   = "fasta"


		qc_input = inFile

		#start by applying filters to R1, if appropriate
		if arguments['--filterOptions'] != "None":
			print("QCing %s" % inFile, file=sys.stderr)
			filter_options = arguments['--filterOptions'].split(" ")
			subprocess.call([vsearch,
					 '-fastx_filter', inFile,
					 '--fastqout', "%s/r1_f%d_filtered.fq"%(prj_tree.preprocess,fileNum)] +
					filter_options, stderr = logFile)

			qc_input = "%s/r1_filtered.fq"%prj_tree.preprocess

			#check for R2
			if len(arguments['--reverse']) > 0:
				print("QCing %s" % arguments['--reverse'][fileNum], file=sys.stderr)
				subprocess.call([vsearch,
						 '-fastx_filter', arguments['--reverse'][fileNum],
						 '--fastqout', "%s/r2_f%d_filtered.fq"%(prj_tree.preprocess,fileNum)] +
						filter_options, stderr = logFile)
				#now merge
				merge_options = arguments['--mergeOptions'].split(" ")
				subprocess.call([vsearch,
						 '-fastq_mergepairs', "%s/r1_f%d_filtered.fq"%(prj_tree.preprocess,fileNum),
						 '-reverse', "%s/r2_f%d_filtered.fq"%(prj_tree.preprocess,fileNum),
						 '--fastqout', "%s/f%d_merged.fq"%(prj_tree.preprocess,fileNum)] +
						merge_options, stderr = logFile)

				qc_input = "%s/f%d_merged.fq"%(prj_tree.preprocess,fileNum)

		elif len(arguments['--reverse']) > 0:
			#merge without QC
			merge_options = arguments['--mergeOptions'].split(" ")
			subprocess.call([vsearch,
					 '-fastq_mergepairs', arguments['--input'][fileNum],
					 '-reverse', arguments['--reverse'][fileNum],
					 '--fastqout', "%s/f%d_merged.fq"%(prj_tree.preprocess,fileNum)] +
					merge_options, stderr = logFile)

			qc_input = "%s/f%d_merged.fq"%(prj_tree.preprocess,fileNum)

		if arguments['--printQC'] is not None:
			print( "Calculating FastQ stats...", file=sys.stderr)
			subprocess.call([vsearch,
					 '-fastq_stats', qc_input,
					 '--log', arguments['--printQC']], #this presumably overwrites if multiple inputs are specified
					stderr = logFile)

		#dereplicate the input to save time
		#want to preserve qual info for UMI processing - will reimplement if/when vsearch adds fastx_uniques command
		#subprocess.call([vsearch, '-derep_fulllength', qc_input, '-sizeout', '-output', "%s/derep.fa"%prj_tree.preprocess], stderr=logFile)

		processedFiles.append(qc_input)

	#now start processing for umis as long as at least one is defined
	if arguments['--cell'] is not None or arguments['--umi'] is not None or arguments['--r2umi'] is not None:

		#split input into managable chunks
		multiInd = 0
		for myFile in processedFiles:
			with open(myFile, 'r') as toSplit:
				splitForID = SeqIO.parse( toSplit, format=fileformat )
				for chunk in iterator_slice(splitForID, 50000):
					multiInd += 1
					with open( "%s/chunk%04d.%s"%(prj_tree.preprocess, multiInd, fileformat), 'w' ) as writeChunk:
						SeqIO.write( chunk, writeChunk, fileformat )

		#construct the command for find_umis
		umiOpts = ""
		for opt in ['--cell', '--umi', '--r2umi', '--cellWhiteList', '--cellPattern', '--umiWhiteList', '--umiPattern', '--umi2WhiteList', '--umi2Pattern', '--minQ' ]:
			if arguments[opt] is not None:
				umiOpts += " %s '%s'" % (opt, arguments[opt])

		#call find_umis either on cluster or locally
		if arguments['--cluster']:
			with open("%s/umijob.sh"%prj_tree.preprocess, 'w') as jobHandle:
				jobHandle.write(f"#!/bin/bash\n#$ -N findUMIs\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/find_umis.py {prj_tree.preprocess}/chunk$NUM.{fileformat} {fileformat} {umiOpts}\n\n")
			subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%multiInd, "%s/umijob.sh"%prj_tree.preprocess])
		else:
			partial_finder = partial( callFinder, format=fileformat, args=umiOpts )
			pool = Pool(arguments['--threads'])
			pool.map( partial_finder, range(1,multiInd+1) )
			pool.close()
			pool.join()

		#collect output of find_umis
		umi_dict = {}
		pickleFiles = glob.glob(f"{prj_tree.preprocess}/chunk*.pickle")
		random.shuffle( pickleFiles ) #try to keep the load balanced when we split up UMIs for clustering below
		maxUmiReads = 0 #try to estimate memory requirements for clustering
		for p in pickleFiles:
			with open(p, 'rb') as pickle_in:
				chunk_dict = pickle.load(pickle_in)
				for (cb, mi) in chunk_dict:
					if (cb, mi) not in umi_dict:
						umi_dict[ (cb, mi) ] = chunk_dict[ (cb, mi) ].copy()
					else:
						umi_dict[ (cb, mi) ]['count'] += chunk_dict[ (cb, mi) ]['count']
						umi_dict[ (cb, mi) ]['seqs']  += chunk_dict[ (cb, mi) ]['seqs']
					if umi_dict[ (cb, mi) ]['count'] > maxUmiReads:
						maxUmiReads = umi_dict[ (cb, mi) ]['count']

		print("Total: %d sequences in %d UMIs" % ( sum([subdict['count'] for subdict in umi_dict.values()]), len(umi_dict)) )

		#if UMIs are present, generate UMI consensus
		if arguments['--umi'] is not None or arguments['--r2umi'] is not None:
			#print out some details that might be useful for QC
			with open("%s/umi_stats.tsv"%prj_tree.logs, 'w') as handle:
				for cb, mi in umi_dict:
					handle.write("%s\t%s\t%s\n"%(cb,mi,umi_dict[(cb,mi)]['count']))

			#generate pickles to pass to consensus algorithm
			uInd = 0
			for chunk in iterator_slice(umi_dict.values(), 5000):
				uInd += 1
				with open( f"{prj_tree.preprocess}/umi_cons_in_{uInd:04}.pickle", 'wb') as pickle_out:
					pickle.dump( chunk, pickle_out )

			#delete umi_dict to save memory
			umi_dict = None

			#try to estimate memory requirements
			memNeeded = 1<<math.ceil( maxUmiReads / 150 ).bit_length()
			if memNeeded < 8:
				memNeeded = 8 #no need to go below default allocation

			#spawn subprocesses
			if arguments['--cluster']:
				with open("%s/umicons.sh"%prj_tree.preprocess, 'w') as jobHandle:
					jobHandle.write(f"#!/bin/bash\n#$ -N clusterUMIs\n#$-l h_vmem={memNeeded}G\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/cluster_umis.py {prj_tree.preprocess}/umi_cons_in_$NUM.pickle {arguments['--minReads']} {prj_tree.preprocess} --isFeature\n\n")
				subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%uInd, "%s/umicons.sh"%prj_tree.preprocess])
			else:
				partial_cons = partial( getUmiConsensus, minSize=arguments['--minReads'], workdir=prj_tree.preprocess, feature=True)

				pool = Pool(arguments['--threads'])
				blob = pool.map( partial_cons, range(1,uInd+1) )
				pool.close()
				pool.join()

			#collect output
			reps  = list()
			cells = defaultdict( dict )
			small = 0
			multi = 0

			for p in glob.glob(f"{prj_tree.preprocess}/umi_cons_out_*.pickle"):
				with open(p, 'rb') as pickle_in:
					chunk_dict = pickle.load(pickle_in)
					small += chunk_dict['small']
					multi += chunk_dict['multi']
					for c in chunk_dict['results']:
						if c in cells:
							cells[c]['count'] += chunk_dict['results'][c]['count']
							cells[c]['seqs']  += chunk_dict['results'][c]['seqs']
						else:
							cells[c].update( chunk_dict['results'][c] )
						reps  += chunk_dict['results'][c]['seqs']

			print(datetime.datetime.now())
			print( "UMIs saved: %d (in %d cells)\nUMIs with fewer than %d reads: %d\nUMIs with multiple clusters:%d\n\n" % (len(reps),len(cells), arguments['--minReads'],small,multi), file=sys.stderr )
			print( "UMIs saved: %d (in %d cells)\nUMIs with fewer than %d reads: %d\nUMIs with multiple clusters:%d\n\n" % (len(reps),len(cells), arguments['--minReads'],small,multi), file=logFile )

			#write output
			with open( arguments['--umiOutput'], "w" ) as handle:
				SeqIO.write( reps, handle, "fasta" )

			#now do cell barcodes, if present
			if arguments['--cell'] is not None:

				#generate pickles to pass to consensus algorithm
				cInd = 0
				for chunk in iterator_slice(cells.values(), 5000):
					cInd += 1
					with open( f"{prj_tree.preprocess}/cell_cons_in_{cInd:04}.pickle", 'wb') as pickle_out:
						pickle.dump( chunk, pickle_out )

				#spawn subprocesses
				if arguments['--cluster']:
					with open("%s/cellcons.sh"%prj_tree.preprocess, 'w') as jobHandle:
						jobHandle.write(f"#!/bin/bash\n#$ -N clusterCells\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/cluster_umis.py {prj_tree.preprocess}/cell_cons_in_$NUM.pickle {arguments['--minUMIs']} {prj_tree.preprocess} --isCell\n\n")
					subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%cInd, "%s/cellcons.sh"%prj_tree.preprocess])
				else:
					partial_cons = partial( getUmiConsensus, minSize=arguments['--minUMIs'], workdir=prj_tree.preprocess, isCell=True)

					pool = Pool(arguments['--threads'])
					blob = pool.map( partial_cons, range(1,cInd+1) )
					pool.close()
					pool.join()

				#collect output
				final_seqs  = list()
				small = 0

				for p in glob.glob(f"{prj_tree.preprocess}/cell_cons_out_*.pickle"):
					with open(p, 'rb') as pickle_in:
						chunk_dict = pickle.load(pickle_in)
						small += chunk_dict['small']
						for c in chunk_dict['results']:
							final_seqs  += chunk_dict['results'][c]['seqs']

				print("%s: %d sequences discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']) , file=sys.stderr)
				print("%s: %d sequences discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']) , file=logFile)
				with open( arguments['--cellOutput'], "w" ) as handle:
					SeqIO.write( final_seqs, handle, "fasta" )

			else:
				#UMIs only, do dereplication and collision removal
				pass

		else:
			#no UMIs present, only cell barcodes
			#in this case, minUMIs effectively replaces minReads

			#generate pickles to pass to consensus algorithm
			cInd = 0
			for chunk in iterator_slice(umi_dict.values(), 5000):
				cInd += 1
				with open( f"{prj_tree.preprocess}/cell_cons_in_{cInd:04}.pickle", 'wb') as pickle_out:
					pickle.dump( chunk, pickle_out )

			#spawn subprocesses
			if arguments['--cluster']:
				with open("%s/cellcons.sh"%prj_tree.preprocess, 'w') as jobHandle:
					jobHandle.write(f"#!/bin/bash\n#$ -N clusterCells\n#$-cwd\nNUM=`printf \"%04d\" $SGE_TASK_ID`\n\nmodule load Biopython/1.73-foss-2016b-Python-3.6.7\n\n{SCRIPT_FOLDER}/annotate/cluster_umis.py {prj_tree.preprocess}/cell_cons_in_$NUM.pickle {arguments['--minUMIs']} {prj_tree.preprocess} --isCell\n\n")
				subprocess.call([qsub, '-sync', 'y', '-t', "1-%d"%cInd, "%s/cellcons.sh"%prj_tree.preprocess])
			else:
				partial_cons = partial( getUmiConsensus, minSize=arguments['--minUMIs'], workdir=prj_tree.preprocess, isCell=True)
				pool = Pool(arguments['--threads'])
				blob = pool.map( partial_cons, range(1,cInd+1) )
				pool.close()
				pool.join()

			#collect output
			final_seqs  = list()
			small = 0

			for p in glob.glob(f"{prj_tree.preprocess}/cell_cons_out_*.pickle"):
				with open(p, 'rb') as pickle_in:
					chunk_dict = pickle.load(pickle_in)
					small += chunk_dict['small']
					for c in chunk_dict['results']:
						final_seqs  += chunk_dict['results'][c]['seqs']

			print("%s: %d sequences discarded because they contained fewer than %d reads..." % (datetime.datetime.now(), small, arguments['--minUMIs']), file=sys.stderr)
			print("%s: %d sequences discarded because they contained fewer than %d reads..." % (datetime.datetime.now(), small, arguments['--minUMIs']), file=logFile)
			with open( arguments['--cellOutput'], "w" ) as handle:
				SeqIO.write( final_seqs, handle, "fasta" )

	else:
		#anything special to do if there are no UMIs/barcodes at all?
		temp = processedFiles
		processedFiles = []
		for f in temp:
			shutil.move( f, re.sub(prj_tree.preprocess,prj_tree.home,f) )
			processedFiles.append( re.sub(prj_tree.preprocess,prj_tree.home,f) )

	# call 1.1 if requested
	if arguments['--runVBlast']:
		cmd = "%s/annotate/1.1-blast_V.py" % SCRIPT_FOLDER

		#ok we'll do one little bit of sanity checking, since these mutually exclusive options both have default values
		if arguments['--cluster']:
			arguments['--threads'] = None

		for opt in [ '--locus', '--lib', '--species', '--npf', '--minl', '--maxl',
           '--jlib', '--dlib', '--clib', '--jmotif', '--nterm', '--file',
					 '--min1', '--min2', '--id', '--maxgaps', '--rearrangements', '--threads']:
			if arguments[opt] is not None:
				cmd += " %s '%s'" % (opt, arguments[opt])
		for flag in ['--derep', '--cluster', '-f', '--runJBlast', '--noD', '--noC',
		             '--runFinalize', '--noclean', '--runClustering', '--runCellStatistics']:
		    if arguments[flag]:
		        cmd += " %s" % flag
		if arguments['--cell'] is not None:
			cmd += " --fasta %s" % arguments['--cellOutput']
		elif arguments['--umi'] is not None or arguments['--r2umi'] is not None:
			cmd += " --fasta %s" % arguments['--umiOutput']
		else:
			for fastaFile in processedFiles:
				cmd += " --fasta %s" % fastaFile

		print( "Calling 1.1 with command line: %s" % cmd )
		os.system( cmd )

	# clean up clustering files
	if not arguments["--keepWorkFiles"]:
		to_clean = glob.glob("%s/*"%prj_tree.preprocess)
		if len(to_clean) > 0:
			print("Cleaning up old files (this may take a while)...",file=sys.stderr)
			for f in to_clean:
				try:
					os.remove(f)
				except IsADirectoryError:
					shutil.rmtree(f)


if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--minQ']     = int( arguments['--minQ'] )
	arguments['--minReads'] = int( arguments['--minReads'] )
	arguments['--minUMIs']	= int( arguments['--minUMIs'] )
	arguments['--threads']  = int( arguments['--threads'] )

	if arguments['--cluster']:
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

	if len(arguments['--input']) != len(arguments['--reverse']) and len(arguments['--reverse']) > 0:
		sys.exit( "The --reverse option must be specified the same number of times as --input!" )

	if len(arguments['--featureLibrary']) > 2:
		sys.exit( "The --featureLibrary option takes a maximum of 2 files.")
	if len(arguments['--featureLibrary']) > 0 and (arguments['--cell'] is None or arguments['--umi'] is None):
		sys.exit( "Cell hashing and feature barcoding require both --cell and --umi to be specificed.")

	if not all([os.path.isfile(x) for x in arguments['--input']]):
		sys.exit( "One or more input files are missing" )

	if len(arguments['--reverse']) > 0 and not all([os.path.isfile(x) for x in arguments['--reverse']]):
		sys.exit( "One or more R2 files are missing" )

	prj_tree = ProjectFolders( os.getcwd() )
	prj_name = fullpath2last_folder(prj_tree.home)

	old_files = glob.glob("%s/*"%prj_tree.preprocess) + glob.glob("byUMI.fa") + glob.glob("byCell.fa")
	if len(old_files) > 0:
		if arguments['-f']:
			print("Cleaning up old files...",file=sys.stderr)
			for f in old_files:
				try:
					os.remove(f)
				except IsADirectoryError:
					shutil.rmtree(f)
		else:
			sys.exit( "Old files exist: Please use the -f flag to force the start of a new analysis" )


	#log command line
	logCmdLine(sys.argv)

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


	#open the logfile
	logFile = open( arguments['--logFile'], "w" )

	main()
