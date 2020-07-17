#!/usr/bin/env python3

"""
2.1-calculate_id-div.py

This script uses MUSCLE or ClustalO to calculate sequence identity between 
      reads and the assigned germline V gene, as well as between the read
      and known antibodies of interest. Outputs two files: _coverage.tab,
      which gives the coverage of the reference sequence by the query and
      _id-div.tab, which gives the *divergence* from germline V and the
      *identity* to other reference sequences. The _id-div.tab file is used
      for making I-D plots with 4.3-plot_identity_divergence.R. Germline V
      identity will also be added to the AIRR rearrangements file.

Usage: 2.1-calculate_id-div.py [ -f input.fa -g germlines.fa -a antibodies.fa -o output -t 1 --align muscle --gap mismatch -d ]

Options:
     -f input.fa        Sequence file to be annotated. In order to calulate
                           identity to germline, the assigned V gene must be
                           present in the fasta def line as eg "v_call=IGHV1-2*02"
                           (the old format of "V_gene=" will also still work).
                           [default: output/sequences/nucleotide/<project>_goodVJ_unique.fa]
     -g germline.fa     File with germline genes used in annotation step.
                           [default: <SONAR>/germDB/IgHKLV_cysTruncated.fa]
     -a antibodies.fa   Fasta file with the sequences of known antibodies that the
                           NGS data should be compared to.
     -o output          Specify directory and file stem for output; "_coverage.tab"
                           and "_id-div.tab" will be appended. If not specified,
                           output will be in the output/tables directory (or the
                           current directory if output/tables doesn't exist) and
                           use the same stem as the input file.
     -t 1               Number of threads to use for alignments. [default: 1]
     --align muscle     Program to use for alignments. Options are 'muscle'
                           and 'clustalo'. [default: muscle]
     --gap mismatch     How to count gaps in the alignment. Options are
                           'mismatch' and 'ignore'. [default: mismatch]
     -d                 A flag to indicate that vsearch should be used to deduplicate
                           the input sequences before running the alignments. Can
                           save significant time for large files which have not
                           already been deduplicated or clustered. All input sequences
                           will still appear in the output file. [default: False]

Created by Zizhang Sheng.
Modified to use VSearch by Chaim A Schramm 2018-07-30.
Ported to Python to handle AIRR-formatted data by CAS 2018-10-17.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National 
                         Institutes of Health, USA. All rights reserved.

"""


import sys, os, re, glob
from docopt import docopt
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from multiprocessing import Pool
import airr


try:
	from SONAR.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/lineage")
	sys.path.append(find_SONAR[0])
	from SONAR.lineage import *



def scoreSeqs( refSeq, querySeq, tempFile):

	#sanity check
	refSeq.seq   = re.sub( "-", "", str(refSeq.seq) )
	querySeq.seq = re.sub( "-", "", str(querySeq.seq) )
	
	with open("%s.fa"%tempFile, "w") as handle:
		handle.write(">%s\n%s\n>%s\n%s\n" % ( refSeq.id, refSeq.seq, querySeq.id, querySeq.seq ))

	align_cline = MuscleCommandline(cmd=muscle, input="%s.fa"%tempFile, out="%s.aln"%tempFile,
									diags=True, maxiters=2, gapopen=-1000.0)
	if arguments['--align'] == "clustalo":
		align_cline = ClustalOmegaCommandline(cmd=clustalo, infile="%s.fa"%tempFile, outfile="%s.aln"%tempFile, force=True)
	try:
		stdout, stderr = align_cline()
	except:
		sys.stderr.write( "Error alinging %s to %s (will skip)\n" % (querySeq.id, refSeq.id, stderr) )
		return "NA","NA"

	alignment = AlignIO.read("%s.aln"%tempFile, "fasta")
	
	#Input order is not maintained, so we need a little
	#   kludge to check which one is the reference sequence.
	refRow = 0
	if alignment[0].id != refSeq.id:
		refRow=1

	#trim terminal gaps
	leftGap = re.match( "-+", str(alignment[refRow].seq) )
	if not leftGap:
		leftGap = re.match( "-+", str(alignment[1-refRow].seq) )
	if leftGap:
		alignment = alignment[:, leftGap.end():]
	rightGap = re.search( "-+$", str(alignment[refRow].seq) )
	if not rightGap:
		rightGap = re.match( "-+", str(alignment[1-refRow].seq) )
	if rightGap:
		alignment = alignment[:, 0:rightGap.start()]

	#check coverage of reference sequence
	coverage = 100 * len( re.sub("-", "", str(alignment[refRow].seq) ) ) / len( refSeq.seq )

	match, gap = 0,0
	for a,b in zip( alignment[refRow].seq, alignment[1-refRow].seq ):
		if b == a:
			match += 1
		elif b == "-" or a == "-":
			gap += 1

	ident = 100 * match / alignment.get_alignment_length()
	if arguments['--gap'] == "ignore":
		ident = 100 * match / ( alignment.get_alignment_length() - gap )

	return coverage, ident



def runAlign( fileName ):

	global germs, mature

	sequences = []
	results	  = dict()

	tempName = re.sub( "\.fa", "_temp", fileName )
	
	print( "Starting work on %s..." % fileName )

	#load sequences and preprocess to align:
	reader = SeqIO.parse(open(fileName, "r"), "fasta")
	for entry in reader:
		results[entry.id] = dict()
		gene = re.search("(v_call|V_gene)=((IG[HKL]V[^*]+)[^,\s]+)",entry.description)
		if gene:
			germline = gene.groups()[1]
			results[entry.id]['vlookup'] = gene.groups()[2]

			if not germline in germs:
				sys.stderr.write( "%s might be misassigned; %s is not in my germline library. Skipping...\n" % (entry.id, germline) )
				results[entry.id]['germline'] = ("NA", "NA")
			else:
				results[entry.id]['germline'] = scoreSeqs( germs[germline], entry, tempName )
		else:
			sys.stderr.write( "Error, can't find germline V for %s...\n" % entry.id )
			results[entry.id]['germline'] = ("NA", "NA")
			results[entry.id]['vlookup'] = "unknown"
		
		for nat in mature:
			results[entry.id][nat] = scoreSeqs( mature[nat], entry, tempName )

	for f in glob.glob( "%s.*" % tempName ):
		os.remove( f )

	return results



def main():
	
	global germs
	germs = dict()
	for entry in SeqIO.parse(open(arguments['-g'], "r"), "fasta"):
		germs[entry.id] = entry

	global mature
	mature = dict()
	if arguments['-a'] is not None:
		for entry in SeqIO.parse(open(arguments['-a'], "r"), "fasta"):
			mature[entry.id] = entry

			
	inputFile = arguments['-f']
	dedup	  = dict()
	if arguments['-d']:
		subprocess.call( [ vsearch, "-derep_fulllength", arguments['-f'],
				   "-output", "temp_dedup.fa",
				   "-uc", "temp.uc", "-notrunclabels" ] )

		inputFile = "temp_dedup.fa"

		#process the uc file
		with open("temp.uc" ,"r") as handle:
			uc = csv.reader(handle, delimiter="\t")
			for row in uc:
				if row[0] == "S":
					dedup[ row[8].split(" ")[0] ] = row[8].split(" ")[0]
				elif row[0] == "H":
					dedup[ row[8].split(" ")[0] ] = row[9].split(" ")[0]


	results	  = dict()
	#If we are multithreading, split input into chunks
	if arguments['-t'] > 1:
		index	= 0
		counter = 0
		chunk	= []
		reader	= SeqIO.parse(open(inputFile, "r"), "fasta")
		for entry in reader:
			chunk.append(entry)
			counter += 1
			if counter == 1000:
				with open("%s/align/align%06d.fa" % (prj_tree.lineage, index), "w") as handle:
					SeqIO.write( chunk, handle, "fasta" )
				index += 1
				counter = 0
				chunk = []
		if counter > 0:
			with open("%s/align/align%06d.fa" % (prj_tree.lineage, index), "w") as handle:
				SeqIO.write( chunk, handle, "fasta" )
			index += 1 #so we can use range properly

		#now create a pool and start the actual work
		filterPool = Pool( arguments['-t'] )
		dataBlob = filterPool.map( runAlign, [ "%s/align/align%06d.fa" % (prj_tree.lineage, i) for i in range(index) ] )
		filterPool.close()
		filterPool.join()

		#Recover results
		for blob in dataBlob:
			results.update( blob )
			
	else:
		#unthreaded, just do the whole thing
		results = runAlign(inputFile)

	
	#get some outputs set up
	outFile	 = os.path.basename( os.path.splitext( arguments['-f'] )[0] )
	if os.path.isdir( prj_tree.tables ):
		outFile = "output/tables/" + outFile
	if arguments['-o'] is not None:
		outFile = arguments['-o']

	nats = sorted( mature.keys() )
	
	covFile	 = open("%s_coverage.tab"%outFile, "w")
	coverage = csv.writer( covFile, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE )
	coverage.writerow( ['sequence_id', 'germ_cov'] + nats )
	
	idFile	 = open("%s_id-div.tab"%outFile, "w")
	iddiv	 = csv.writer( idFile, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE )
	iddiv.writerow( ['sequence_id', 'v_gene', 'germ_div'] + nats )

	#sort the freaking list and output
	if arguments['-d']:
		for s in sorted( dedup.keys() ):
			(germc, germi) = results[dedup[s]]['germline']
			if not germc == "NA":
				germc = "%.1f" % germc
				germi = "%.1f" % (100-germi)
			coverage.writerow( [ s, germc ] +
					   [ "NA" if results[dedup[s]][n][0]=="NA" else "%.1f"%results[dedup[s]][n][0] for n in nats ] )
			iddiv.writerow( [ s, results[dedup[s]]['vlookup'], germi ] +
					[ "NA" if results[dedup[s]][n][1]=="NA" else "%.1f"%results[dedup[s]][n][1] for n in nats ] )

		#take this opportunity to do some cleanup
		os.remove("temp_dedup.fa")
		os.remove("temp.uc")
	else:
		for s in sorted( results.keys() ):
			(germc, germi) = results[s]['germline']
			if not germc == "NA":
				germc = "%.1f" % germc
				germi = "%.1f" % (100-germi)
			coverage.writerow( [ s, germc ] +
					   [ "NA" if results[s][n][0]=="NA" else "%.1f"%results[s][n][0] for n in nats ] )
			iddiv.writerow( [ s, results[s]['vlookup'], germi ] +
					[ "NA" if results[s][n][1]=="NA" else "%.1f"%results[s][n][1] for n in nats ] )
	
	covFile.close()
	idFile.close()

	
	#do AIRR output
	if os.path.dirname(arguments['-f']) == "output/sequences/nucleotide" and not 'CDR3' in arguments['-f']:
		if os.path.isfile("%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name)):
			withDiv = airr.derive_rearrangement( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name), fields=['v_identity'] )
			for r in airr.read_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) ):
				if dedup.get(r['sequence_id'], r['sequence_id']) in results:
					# omit NAs here to comply with AIRR format
					if not results[ dedup.get(r['sequence_id'], r['sequence_id']) ]['germline'][1] == "NA":
						r['v_identity'] = "%0.3f" % (results[ dedup.get(r['sequence_id'], r['sequence_id']) ]['germline'][1]/100)
				withDiv.write(r)
			withDiv.close()
			os.rename( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) )

                        

if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['-t'] = int( arguments['-t'] )

	if arguments['--align'] not in ['muscle', 'clustalo']:
		sys.exit( "Error: recognized alignment programs are 'muscle' and 'clustalo' only" )    
	
	if arguments['--gap'] not in ['mismatch', 'ignore']:
		sys.exit( "Error: recognized gap options are 'mismatch' and 'ignore' only" )	
	
	arguments['-g'] = re.sub( "<SONAR>", SCRIPT_FOLDER, arguments['-g'] )

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	arguments['-f'] = re.sub("<project>", prj_name, arguments['-f'])
	if not os.path.isfile(arguments['-f']):
		sys.exit( "Can't find input file %s" % arguments['-f'] )

	#make dir for alignments but ignore error if it already exists
	os.makedirs("%s/align"%prj_tree.lineage, exist_ok=True)
	       
	#log command line
	logCmdLine(sys.argv)	

	main()
