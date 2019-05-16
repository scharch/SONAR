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
   * Come up with a better way to decide which cells are "real"; 1 UMI might be okay if it has 50
         reads; OTOH 3 UMIs might not be enough if they do not all cluster together.

Usage: 1.0-preprocess.py --input read1.fq... [ --reverse read2.fq... ] [ --cellWhiteList barcodes.txt | --cellPattern NNNNNN ] [ --umiWhiteList barcodes.txt | --umiPattern NNNNNN ] [ --umi2WhiteList barcodes.txt | --umi2Pattern NNNNNN ] [ options ]

Options:
    --input read1.fq               File with raw data to process. Can be used multiple times **IF**
                                       files represent technical or sequencing replicates, where
                                       overlapping cell barcodes/UMIs can be assumed to be derived
                                       from the same physical source.
    --reverse read2.fq             File with reverse reads for paired ends. If used, must be provided
                                       the same number of times (and in the same order) as --input.
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
    --minQ                         Minimum PHRED score for all bases in a UMI or cell barcode. Reads 
                                       with *any* base in the UMI/barcode below this threshold will
                                       be discarded. [default: 20]
    --minReads 1                   Minimum number of reads per UMI. UMIs with fewer reads will be
                                       discarded. [default: 3]
    --minUMIs 1                    Minimum number of UMIs per cell. Cells with fewer UMIs will be
                                       discarded. [default: 3]
    --umiOutput file.fa            File in which to save the UMI-processed sequences. Ignored if
                                       UMIs are not present. [default: byUMI.fa]
    --cellOutput file.fa           File in which to save the cell-processed sequences. Ignored if
                                       cell barcodes are not present. [default: byCell.fa]
    --runVBlast                    Flag to call 1.1 script when finished. Additional options to that
                                       script (including the ability to automatically call 1.2 and
                                       further scripts in Module 1) are listed below. This script
                                       will not check the validity of options passed downstream, so
                                       user beware.
    --threads 1                    Number of threads to use. [default: 1]
    -f                             Flag to force overwriting of old files. [default: False]

Options for other annotation scripts (see those help messages for details):
    --locus H
    --lib LIB
    --derep
    --minl <300>
    --maxl <600>
    --cluster
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

Copyright (c) 2019 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os, shutil, gzip
from docopt import docopt
import itertools
from collections import defaultdict
from Bio import SeqIO
from multiprocessing import Pool
import datetime
from functools import partial

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *
	

# a utility function to get us a slice of an iterator, as an iterator
# when working with iterators maximum lazyness is preferred
# from https://stackoverflow.com/a/44502827
def iterator_slice(iterator, length):
	iterator = iter(iterator)
	slice = 0
	while True:
		res = tuple(itertools.islice(iterator, length))
		if not res:
			break
		slice += 1
		yield slice, res


def getUmiConsensus(iter_tuple, minSize, workdir, isUMI=True):

	chunk, umi_iter = iter_tuple

	results = {}
	small	= 0
	multi	= 0
	
	print( "Generating consensus on chunk #%d..." % chunk)

	for umi in umi_iter:

		#check read threshold
		if umi['count'] < minSize:
			small += 1
			continue

		if len(umi['seqs']) == 1:
			#save time on singletons (if they weren't excluded by the read threshold)
			if isUMI:
				umi['seqs'][0].id += ";seqs=1;size=%d;consensus_count=%d" % (umi['count'],umi['count'])
				umi['seqs'][0].description = ""
			else:
				umi['seqs'][0].id = "%s.1 cell_id=%s duplicate_count=1 consensus_count=%s"%( umi['cell'], umi['cell'], umi['count'] )
				umi['seqs'][0].description = ""

			if (umi['cell']) not in results:
				results[ umi['cell'] ] = { 'cell':umi['cell'], 'umi':umi['cell'], 'count':1, 'seqs':umi['seqs'].copy() }
			else:
				results[ umi['cell'] ]['count'] += 1
				results[ umi['cell'] ]['seqs']  += umi['seqs']

		else:

			subdir = workdir
			#iddef  = "3"
			if isUMI:
				subdir += "/%s" % umi['cell']
				#iddef = "2"

			#cluster and rapid align with vsearch
			os.makedirs(subdir, exist_ok=True)
			with open("%s/%s.fa" % (subdir, umi['umi']), "w") as handle:
				SeqIO.write(umi['seqs'], handle, "fasta")

			subprocess.call([vsearch,
					 "-cluster_fast", "%s/%s.fa" % (subdir, umi['umi']),
					 "-consout", "%s/%s_cons.fa" % (subdir, umi['umi']),
					 "-id", "0.95",
					 "-iddef", "3",
					 "-sizein", "-sizeout",
					 "-mincols", '150',
					 "-gapopen", "10I/10E", #lower gap open penalty to better account for internal indels
					 "-gapext", "2I/2E", #don't make endgaps cheaper; encourages TSOs to align properly
					 "-clusterout_sort", #so we can look at just the biggest
					 "-quiet" #supress screen clutter
					 ])

			with open("%s/%s_cons.fa" % (subdir, umi['umi']), 'r') as cons_file:
				seq_number = 0
				for cons in SeqIO.parse(cons_file, "fasta"):
					seq_number += 1
					if isUMI and seq_number > 1:
						#how to handle more than one cluster per umi?
						#  -depends on presence/absence of cell barcodes, I guess. user param?
						#use it to do error checking???
						multi += 1
						break

					num_reads = re.search(";seqs=(\d+);size=(\d+)",cons.id)
					if num_reads:
						if isUMI:
							if int(num_reads.group(2)) < minSize:
								small += 1
								continue
							else:
								cons.id += ";consensus_count=%s" % num_reads.group(2) #save size annotation for further clustering/dereplication
						else:
							cons.id	 = "%s.%d cell_id=%s duplicate_count=%s consensus_count=%s"%( umi['cell'], seq_number, umi['cell'], num_reads.group(1), num_reads.group(2) )				
						
						cons.description = ""
						if (umi['cell']) not in results:
							results[ umi['cell'] ] = { 'cell':umi['cell'], 'umi':umi['cell'], 'count':1, 'seqs':[cons] }
						else:
							results[ umi['cell'] ]['count'] += 1
							results[ umi['cell'] ]['seqs'].append(cons)

	return (small, multi, results)


def main():

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
		total_count = 0

		for myFile in processedFiles:

			file_umis = {}
			count     = 0
			bad_umi   = 0
			low_qual  = 0
			print("%s: Starting to look for UMIs in %s" % (datetime.datetime.now(), myFile), file=sys.stderr)

			with open(myFile, "r") as handle:
				for seq in SeqIO.parse( handle, fileformat):
					count += 1
					if count % 50000 == 0: print( "%s: Processed %d sequences in %d UMIs so far; Discarded %d reads with low quality UMIs and %d additional reads with illegal UMIs." % (myFile, count, len(file_umis), low_qual, bad_umi) )

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
						if fileformat=="fastq" and any([ x<arguments['--minQ'] for x in seq.letter_annotations['phred_quality'][cb_start:cb_end] ]):
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
						if fileformat=="fastq" and any([ x<arguments['--minQ'] for x in seq.letter_annotations['phred_quality'][umi_start:umi_end] ]):
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
						if fileformat=="fastq" and any([ x<arguments['--minQ'] for x in seq.reverse_complement().letter_annotations['phred_quality'][umi2_start:umi2_end] ]):
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

					if cell_barcode != "":
						seq.id += ";cell=%s"%cell_barcode
					if molecule_id	!= "":
						seq.id += ";umi=%s"%molecule_id

					if (cell_barcode, molecule_id) not in umi_dict:
						umi_dict[ (cell_barcode, molecule_id) ] = { 'cell':cell_barcode, 'umi':molecule_id, 'count':reads, 'seqs':[seq] }
					else:
						umi_dict[ (cell_barcode, molecule_id) ]['count'] += reads
						umi_dict[ (cell_barcode, molecule_id) ]['seqs'].append(seq)

					if (cell_barcode, molecule_id) not in file_umis:
						file_umis[ (cell_barcode, molecule_id) ] = 1
						
					#next sequence in file
					
			print( "Finished %s: %d sequences in %d UMIs; Discarded %d reads with low quality UMIs and %d additional reads with illegal UMIs." % (myFile, count, len(file_umis), low_qual, bad_umi) )
			total_count += count
			#loop to next input file, if relevant

			
		print("Total: %d sequences in %d UMIs" % (total_count, len(umi_dict)) )
		      
		#if UMIs are present, generate UMI consensus
		if arguments['--umi'] is not None or arguments['--r2umi'] is not None:
			#print out some details that might be useful for QC
			with open("%s/umi_stats.tsv"%prj_tree.logs, 'w') as handle:
				for cb, mi in umi_dict:
					handle.write("%s\t%s\t%s\n"%(cb,mi,umi_dict[(cb,mi)]['count']))
					
			#go through each lineage and do the alignment
			reps  = list()
			cells = defaultdict( dict )
			small = 0
			multi = 0

			partial_cons = partial( getUmiConsensus, minSize=arguments['--minReads'], workdir=prj_tree.preprocess, isUMI=True)

			pool = Pool(arguments['--threads'])
			blob = pool.map( partial_cons, iterator_slice(umi_dict.values(), 1000) ) #number per slice needs optimization
			pool.close()
			pool.join()

			for s,m,r in blob:
				small += s
				multi += m
				for c in r:
					if c in cells:
						cells[c]['count'] += r[c]['count']
						cells[c]['seqs']  += r[c]['seqs']
					else:
						cells[c].update( r[c] )
					reps  += r[c]['seqs']

			print(datetime.datetime.now())
			print( "UMIs saved: %d (in %d cells)\nUMIs with fewer than %d reads: %d\nUMIs with multiple clusters:%d\n\n" % (len(reps),len(cells), arguments['--minReads'],small,multi), file=sys.stderr )
			print( "UMIs saved: %d (in %d cells)\nUMIs with fewer than %d reads: %d\nUMIs with multiple clusters:%d\n\n" % (len(reps),len(cells), arguments['--minReads'],small,multi), file=logFile )
	
			#write output (make this a user param)
			with open( arguments['--umiOutput'], "w" ) as handle:
				SeqIO.write( reps, handle, "fasta" )

			#now do cell barcodes, if present
			if arguments['--cell'] is not None:
				final_seqs  = list()
				small = 0
				
				partial_cons = partial( getUmiConsensus, minSize=arguments['--minUMIs'], workdir=prj_tree.preprocess, isUMI=False)

				pool = Pool(arguments['--threads'])
				blob = pool.map( partial_cons, iterator_slice(cells.values(), 1000) ) #number per slice needs optimization
				pool.close()
				pool.join()

				for s,m,r in blob:
					small += s
					for c in r:
						final_seqs  += r[c]['seqs']

				print("%s: %d cells discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']) , file=sys.stderr)
				print("%s: %d cells discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']) , file=logFile)
				with open( arguments['--cellOutput'], "w" ) as handle:
					SeqIO.write( final_seqs, handle, "fasta" )

			else:
				#UMIs only, do dereplication and collision removal
				pass

		else:
			#no UMIs present, only cell barcodes
			#in this case, minUMIs effectively becomes a minReads per cell,
			#    with no minimum on the reads per final sequence.
			#Not clear this is the desired behavior, may need to add a
			#    third case to the getUmiConsensus function.
			final_seqs  = list()
			small = 0
				
			partial_cons = partial( getUmiConsensus, minSize=arguments['--minUMIs'], workdir=prj_tree.preprocess, isUMI=False)

			pool = Pool(arguments['--threads'])
			blob = pool.map( partial_cons, iterator_slice(umi_dict.values(), 1000) ) #number per slice needs optimization
			pool.close()
			pool.join()

			for s,m,r in blob:
				small += s
				for c in r:
					final_seqs  += r[c]['seqs']

			print("%s: %d cells discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']), file=sys.stderr)
			print("%s: %d cells discarded because they contained fewer than %d UMIs..." % (datetime.datetime.now(), small, arguments['--minUMIs']), file=logFile)
			with open( arguments['--cellOutput'], "w" ) as handle:
				SeqIO.write( final_seqs, handle, "fasta" )

	else:
		#anything special to do if there are no UMIs/barcodes at all?
		pass

	# call 1.1 if requested
	if arguments['--runVBlast']:
		cmd = "%s/annotate/1.1-blast_V.py" % SCRIPT_FOLDER

		#ok we'll do one little bit of sanity checking, since this script doesn't use the --cluster flag
		if arguments['--cluster']:
			arguments['--threads'] = None
		
		for opt in [ '--locus', '--lib', '--npf', '--minl', '--maxl', '--jlib', '--dlib', 
					 '--clib', '--jmotif', '--nterm', '--file', '--min1', '--min2',
					 '--id', '--maxgaps', '--rearrangements', '--threads']: 
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
	
	if len(arguments['--input']) != len(arguments['--reverse']) and len(arguments['--reverse']) > 0:
		sys.exit( "The --reverse option must be specified the same number of times as --input!" )
		
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

