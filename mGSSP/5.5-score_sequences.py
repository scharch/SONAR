#!/usr/bin/env python3

'''
5.5-score_sequences.py

This script scores the rarity for each position in a set of antibody sequences
    and identify those qualifying as "extremely rare."

Usage: 5.5-score_sequences.py QVQLVQ... -v IGHV1-2 [ --gssp GSSP.txt ] [options]
       5.5-score_sequences.py -f input.fasta [ -v IGHV1-2 --gssp GSSP.txt ] [options]
       5.5-score_sequences.py -r rearrangements.tsv [ --gssp GSSP.txt ] [options]

Options:
    QVQLVQ                   Amino acid sequence of one or more antibodies to score.
    -f input.fasta           Fasta file containing sequences to score.
    -r rearrangements.tsv    AIRR-formatted rearrangments with sequences to score. Assumes -n.
    -v IGHV1-2               V gene of the input sequences. If not provided, will be extracted
                                 from the 'V_gene=' or 'v_call=' tag in the fasta def line or
                                 the 'v_call' field of a rearrangements file.
    --gssp GSSP.txt          File with GSSPs to use for scoring the sequences.
                                 [default: <SONAR>/sample_data/GSSPs/Sheng2017_VH_GSSPs.txt]
    -n                       Flag to indicate that input sequences are nucleotide and must
                                 be translated. [default: False]
    --germ germline.fa       File with germline sequences, used for aligning input before
                                 scoring. [default: <SONAR>/germDB/IgHKLV_cysTruncated.AA.fa]
    --rare .995              Threshold for what counts as 'extremely rare'. [default: 0.995]
    --lineage                Instead of reporting rare substitutions in each sequence, report
    						     those that appear in at least `--threshold` percent of all
    						     input sequences. [default: False]
    --threshold 50           If `--lineage` is specified, the percent of all input sequences
                                 that must contain a rare substitution for it to be reported.
                                 [default: 50]


Created by Chaim A Schramm on 2019-02-22.

Copyright (c) 2019, Vaccine Research Center, National Institutes of Health, USA.
                         All rights reserved.

'''

import sys, csv, re
from docopt import docopt
from Bio import Seq, SeqIO
import airr
from collections import defaultdict

try:
	from SONAR.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/mGSSP")
	sys.path.append(find_SONAR[0])
	from SONAR.mGGSP import *


def checkGermSeq( gene, lib ):
	try:
		germSeq = lib[gene+"*01"]
		return True
	except KeyError:
		print( "Can't find sequence of %s*01 in database %s" % (gene,arguments['--germ']), file=sys.stderr )
		return False
	

def checkGSSP( gene, lib ):
	if not gene in lib:
		print( "Can't find a GSSP for %s in database %s" % (gene,arguments['--gssp']), file=sys.stderr )
		return False
	else:
		return True


def score( sequence, germline, v_rarity ):
	if arguments['-n']:
		try:
			sequence = re.sub("-","",sequence)
		except TypeError:
			#SeqRecord from fasta input
			sequence = re.sub("-","",str(sequence.seq))
		sequence = str(Seq.Seq(sequence).translate())

	align = quickAlign( germline, sequence )
	rare = []
	ind=0
	for r,t in zip(align['ref'],align['test']):
		if r=="-":
			continue
		if ind in v_rarity and t in v_rarity[ind]['mutants'] and v_rarity[ind]['mutants'][t]['average'] >= arguments['--rare']:
			rare.append( f"{r}{ind+1}{t}" )
		ind += 1

	return rare


def main():

	#start by reading in the GSSP
	gssp = GSSP( arguments['--gssp'] )
	gssp.computeRarity()

	#now get germline genes
	germDB = load_fastas( arguments['--germ'] )

	rareSubs = dict()

	if arguments['-r'] is not None:
		for seq in airr.read_rearrangement( arguments['-r'] ):
			gl = re.sub("\*.*","",seq['v_call'])
			if checkGermSeq(gl, germDB) and checkGSSP(gl, gssp.rarity):
				rareSubs[ seq['sequence_id'] ] = score( seq['sequence_alignment'], germDB[gl+"*01"], gssp.rarity[gl] )
	elif arguments['-f'] is not None:
		for seq in generate_read_fasta( arguments['-f'] ):
			if arguments['-v'] is not None:
				if checkGermSeq(arguments['-v'], germDB) and checkGSSP(arguments['-v'], gssp.rarity):
					rareSubs[ seq.id ] = score( seq, germDB[arguments['-v']+"*01"], gssp.rarity[arguments['-v']] )
			else:
				gl = re.search("(v_call|V_gene)=([^\*\s]+)", seq.description)
				if gl:
					if checkGermSeq(gl.group(2), germDB) and checkGSSP(gl.group(2), gssp.rarity):
						rareSubs[ seq.id ] = score( seq, germDB[gl.group(2)+"*01"], gssp.rarity[gl.group(2)] )
				else:
					print("Could not find V gene annotation for %s, skipping..." % seq.id, file=sys.stderr)
					continue
	else:
		if checkGermSeq(arguments['-v'], germDB) and checkGSSP(arguments['-v'], gssp.rarity):
			for sequence in arguments['QVQLVQ']:
				rareSubs[ sequence ] = score( sequence, germDB[arguments['-v']+"*01"], gssp.rarity[arguments['-v']] )
		else:
			sys.exit(1)

	#now do output
	count = 0
	if arguments['--lineage']:
		reverse_dict = defaultdict(list)
		for seq in rareSubs:
			for sub in rareSubs[seq]:
				reverse_dict[ sub ].append( seq )
		for sub in sorted(reverse_dict.keys(), key=lambda x: int(re.search("(\d+)",x).group(1))):
			if 100*len(reverse_dict[sub])/len(rareSubs) >= arguments['--threshold']:
				print(sub)
				count +=1
	else:
		for seq in rareSubs:
			if len(rareSubs[seq]) > 0:
				print( seq + ": " + ",".join(rareSubs[seq]) )
				count +=1

	if count == 0:
		print( "No rare substitutions were found")
		

if __name__ == "__main__":

	arguments = docopt(__doc__)
	
	arguments['--gssp']      = re.sub( "<SONAR>", SCRIPT_FOLDER, arguments['--gssp'] )
	arguments['--germ']      = re.sub( "<SONAR>", SCRIPT_FOLDER, arguments['--germ'] )
	arguments['--rare']	     = float( arguments['--rare'] )
	arguments['--threshold'] = float( arguments['--threshold'] )

	if arguments['-r'] is not None:
		arguments['-n'] = True
	
	#log command line
	logCmdLine(sys.argv)	

	main()
