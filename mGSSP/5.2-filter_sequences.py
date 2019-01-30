#!/usr/bin/env python3

"""

5.2-filter_sequences.py

This script uses pairwise clustal alignments to the assigned germline to detect
      frameshifts in  BCR sequences for use in GSSPs. It will either discard
      these sequences as presumed sequencing error (for functional repertoires)
      or introduce gaps to correct the frameshift (nonfunctional repertoires).
      The script will then translate the sequences and check for at least 
      1 aa substitution from germline; these will be output for building GSSPs.
      Uses clustalw despite it being slower and less efficient because I
      was unable to optimize parameters for muscle or clustalo.

Usage: 5.2-filter_sequences.py IN OUT [ -t 1 --db germ.fa --keep ]

Options
    IN             Fasta file containing sequences to be checked. V gene assignment
                       must be present in def line in standard format.
    OUT            Fasta file in which to save good sequences.
    -t 1           Number of threads to use. [default: 1]
    --db germ.fa   Fasta file containing germline V gene sequences to align to.
                       Truncation of 3' end of V gene is recommended to avoid the
                       introduction of spurious gaps into the alignment due to SHM
                       in CDR3. [default: <sonar>/germDB/IgHKLV_cysTruncated.fa]
    --keep         Flag to indicate that sequences with apparent frameshifts should
                       be retained (ie for passenger repertoires). [default: False]

Created as checkForFrameShift.py by Chaim Schramm on 2015-12-29.
Adapted to align nonfunctional sequences by CAS 2018-07-13.
Added parallel handling of functional sequences, renamed again, and added to mGSSP
        workflow by CAS 2018-08-31.
Multithreaded by CAS 2018-09-05.

Copyright (c) 2015-2018 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""

import glob, os, re, sys
from docopt import docopt
from multiprocessing import Pool
from Bio import Seq
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Alphabet
import numpy

try:
	from sonar.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/mGSSP")
	sys.path.append(find_SONAR[0])
	from sonar.mGSSP import *


def runClustal( fileName ):

	global germs

	sequences   = []
	total, good, mutated = 0, 0, 0

	tempName = re.sub( "\.fa", "_temp", fileName )
	
	print( "Starting work on %s..." % fileName )

	#load sequences and preprocess to align:
	reader = SeqIO.parse(open(fileName, "r"), "fasta")
	for entry in reader:

		total += 1
		gene = re.search("(?:v_call|V_gene)=(IG[HKL]V[^,\s]+)",entry.description)

		if gene:

			germline = gene.groups()[0]

			if not germline in germs:
				print( "%s might be misassigned; %s is not in my germline library. Skipping..." % (entry.id, germline) )
				continue

			with open("%s.fa"%tempName, "w") as handle:
				handle.write(">%s\n%s\n>%s\n%s\n" % ( germline, germs[germline], entry.id, entry.seq ))

			clustal_cline = ClustalwCommandline(cmd=clustalw, infile="%s.fa"%tempName)
			try:
				stdout, stderr = clustal_cline()
			except:
				print( "Error in alignment of %s (will skip): %s" % (entry.id, stderr) )
				for f in glob.glob("%s.*"%tempName): os.remove( os.path.abspath(f) )
				continue

			alignment = AlignIO.read("%s.aln"%tempName, "clustal")

			if not arguments['--keep']:
				shift = False
				for record in alignment:
					#strip end gaps, they don't matter
					#full-codon indels are also fine
					codons = re.sub( "---", "", str(record.seq.strip("-")) )
					if "-" in codons:
						shift = True #likely frameshift --discard!
				if shift:
					for f in glob.glob("%s.*"%tempName): os.remove( os.path.abspath(f) )
					continue
				
			#count no-frameshift seqs for functional repertoires
			good += 1
			
			#now, remove gaps
			
			#Input order is not maintained, so we need a little
			#   kludge to check which one isthe germline sequence.
			germRow = 0
			if alignment[0].id != germline:
				germRow=1

			#look for gaps one at a time so we don't get tripped up by shifting indices
			gap = re.search( "-+", str(alignment[germRow].seq) )
			while (gap):
				alignment = alignment[:, 0:gap.start()] + alignment[:, gap.end():]
				gap = re.search( "-+", str(alignment[germRow].seq) )
					
			#translate to check for AA substitutions
			mySeq = alignment[ 1 - germRow ]
			mySeq.seq, n = re.subn( "-", "N", str(mySeq.seq) )
			mySeq.seq = Seq.Seq( mySeq.seq, Alphabet.IUPAC.ambiguous_dna )
			mySeq.seq = mySeq.seq.translate()
			germSeq = str( alignment[germRow].seq.translate() )

			mutCount = 0
			for a,b in zip( germSeq, str(mySeq.seq) ):
				if b != 'X' and b != "-" and b != a:
					mutCount += 1

			if mutCount > 0:
				mutated += 1
				entry.seq = Seq.Seq( re.sub( "^X+", "", str(mySeq.seq) ), Alphabet.IUPAC.ExtendedIUPACProtein )
				sequences.append(entry)
		
			for f in glob.glob("%s.*"%tempName):
				os.remove( os.path.abspath(f) )

	return (total, good, mutated, sequences)



def main():

	global germs
	germs	    = dict()
	for entry in SeqIO.parse(open(arguments['--db'], "r"), "fasta"):
		germs[entry.id] = entry.seq


	#If we are multithreading, split input into chunks
	if arguments['-t'] > 1:
		index	= 0
		counter = 0
		chunk	= []
		reader	= SeqIO.parse(open(arguments['IN'], "r"), "fasta")
		for entry in reader:
			chunk.append(entry)
			counter += 1
			if counter == 1000:
				with open("%s/work/mGSSP/filter%06d.fa" % (prj_tree.home, index), "w") as handle:
					SeqIO.write( chunk, handle, "fasta" )
				index += 1
				counter = 0
				chunk = []
		if counter > 0:
			with open("%s/work/mGSSP/filter%06d.fa" % (prj_tree.home, index), "w") as handle:
				SeqIO.write( chunk, handle, "fasta" )
			index += 1 #so we can use range properly

		#now create a pool and start the actual work
		filterPool = Pool( arguments['-t'] )
		dataBlob = filterPool.map( runClustal, [ "%s/work/mGSSP/filter%06d.fa" % (prj_tree.home, i) for i in range(index) ] )
		filterPool.close()
		filterPool.join()

		#there must be a better way to recover this data, but using multiprocessing.Value looks complicated
		total,good,mutated=0,0,0
		sequences = []
		for blob in dataBlob:
			total	  += blob[0]
			good	  += blob[1]
			mutated	  += blob[2]
			sequences += blob[3]

	else:
		#unthreaded, just do the whole thing
		total, good, mutated, sequences = runClustal(arguments['IN'])


	SeqIO.write(sequences, arguments['OUT'], "fasta")
	print( "Total: %d, Good: %d, Mutated: %d" % (total, good, mutated) )


if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--db'] = re.sub( "<sonar>", SCRIPT_FOLDER, arguments['--db'] )
	arguments['-t']	  = int( arguments['-t'] )

	#log command line
	logCmdLine(sys.argv)	

	prj_tree = ProjectFolders(os.getcwd())
	os.makedirs("%s/work/mGSSP"%prj_tree.home, exist_ok=True)

	main()
