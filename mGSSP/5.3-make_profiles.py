#!/usr/bin/env python3

"""

5.3-make_profiles.py

This script takes input sequences and generates germline mutability profiles
      for downstream analysis. The recommended workflow is to run the annotation
      module (1.1 through 1.4) followed by 2.4-cluster_into_groups.py and (optionally)
      5.1-repick_lineage_representatives.py. Then filter to remove reads with frameshift
      errors and those that have no amino acid substitutions using 5.2-filter_sequences.py.
      Then run this program.

Usage: 5.3-make_profiles.py <sequences.fa> [ -o profiles.txt -n 300 -p 0 -m 0 -g germV.fa -t 1 -a ]

Options:
   -h --help                    Show this documentation
   <sequences.fa>               Processed ngs sequences to be used to build the profiles
   -o --output profiles.txt     Where to save output [default: profiles.txt]
   -n --numSequences 300        Number of reads to use in building each profile [default: 300]
   -p --profiles 0              Number of profiles to build for each germline gene by randomly 
                                   subsetting -n sequences each time. Currently does not check 
                                   total number of sequences to make sure subsets are different
                                   enough from each other. If set to zero, all sequences will
                                   be used for a single profile, with -n sequences being a
                                   minimum only. [default: 0]
   -m --mask 0                  Mask the first m positions at the 5' end of each gene (missing
                                   data and/or degenerate primers). Can also be specified as a
                                   tsv file with the first column as the gene name and the second
                                   column as the number of positions to mask. In this case, genes
                                   that do not appear in the file will not be masked. [default: 0]
   -g --germline germV.fa       Location of germline V sequences to use for building profiles. 
                                   Expected as trimmed/padded AA sequences. [default: SONAR/germDB/IgHKLV_cysTruncated.AA.fa]
   -t 1                         Number of threads to use. [default: 1]
   -a                           Input sequences are amino acid (don't translate) [default: False]

Created by Chaim Schramm on 2016-05-27.
Added to SONAR as part of mGSSP on 2017-02-24.
Changed some options and defaults by CAS 2018-07-10.
Edited to use Py3 by CAS 2018-08-29.
Renamed to 5.3 and multithreaded by CAS 2018-09-05.
Tweaks for AIRR-formats naming conventions by CAS 2018-10-18.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys, numpy, re, csv, glob, os
from multiprocessing import Pool
from collections import defaultdict, Counter
from docopt import docopt
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

try:
	from SONAR.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/mGSSP")
	sys.path.append(find_SONAR[0])
	from SONAR.annotate import *


def buildGSSP( vgene ):

	results = []

	if len(masterList[vgene]) < arguments["--numSequences"]:
		print( "Skipping %s, not enough sequences (%d)..." % ( vgene, len(masterList[vgene]) ) )
		return []
		
	if vgene not in germList:
		print( "Skipping %s, it's not in the germline database..." %vgene )
		return []

	# Take random overlapping subsets to generate multiple profiles
	#  need to add back a sanity check for capping the number of subsets if there's not enough raw data.
	numProfiles = arguments['--profiles']
	if arguments["--profiles"] == 0:
		numProfiles = 1

	success = 0
		
	for i in range(numProfiles):
		seqs = [] + germList[vgene] #force a copy rather than an alias
		if arguments["--profiles"] == 0:
			seqs += list(masterList[vgene])
		else:
			#get our sequence subset, add the germlines, and write them
			#   to a temporary file for alignment
			seqs += list(numpy.random.choice(masterList[vgene], size=arguments["--numSequences"], replace=False))

		tempFile = "%s/work/mGSSP/%s_profileBuilder" % (prj_tree.home, vgene)
		with open("%s.fa"%tempFile, "w") as temp:
			SeqIO.write(seqs,temp,"fasta")

		muscle_cline = MuscleCommandline(cmd=muscle, input="%s.fa"%tempFile, out="%s.aln"%tempFile)

		#try to speed up the process a little bit for large datasets
		#still going to max out at ~50k seqs per profile (probably)
		muscle_cline.maxiters	= 2
		muscle_cline.diags	= True

		try:
			stdout, stderr = muscle_cline()
		except:
			print( "Error in alignment #%d for %s (skipping)" % (i+1, vgene) )
			for f in glob.glob("%s.*"%tempFile): 
				os.remove(f)
			continue

		alignment = AlignIO.read("%s.aln"%tempFile, "fasta")#"clustal")
		success += 1

		#Input order is not maintained, so we need a little
		#   kludge to find a germline sequences. Use the 
		#   first one to remove any insertions from the alignment
		germRow = 0
		for n, rec in enumerate(alignment):
			if rec.id in [g.id for g in germList[vgene]]:
				germRow = n
				break

		#look for gaps one at a time so we don't get tripped up by shifting indices
		gap = re.search( "-+", str(alignment[germRow].seq) )
		while (gap):
			alignment = alignment[:, 0:gap.start()] + alignment[:, gap.end():]
			gap = re.search( "-+", str(alignment[germRow].seq) )
		
		#Now we get BioPython to make a PSSM for us. To convert that into
		#    a mutability profile, we will delete the germline residue[s]
		#    at each position (but save what they were)
		germRes = defaultdict(Counter)
		summary_align = AlignInfo.SummaryInfo(alignment)
		pssm = summary_align.pos_specific_score_matrix(chars_to_ignore=['-','X'])

		#get number of datapoints at each position (might be different than the number of sequences in the profile if there are gaps or missing data
		# do this by using sum(pos.values()) after ignoring missing data (previous line) but before dumping germline residues.
		denominator = []
		for p,pos in enumerate(pssm):
			denominator.append( sum(pos.values()) - len(germList[vgene]) )
    
		for germ in germList[vgene]:
			for pos, residue in enumerate(germ):
				if residue == "X":
					continue
				germRes[pos][residue] += 1
				pssm[pos][residue] = 0

		#normalize and save
		for p, pos in enumerate(pssm):
			germAA = ",".join([ x[0] for x in germRes[p].most_common() ])
			results.append( [ vgene, i+1, p+1, germAA, "None" if (p < mask[vgene] or denominator[p] < arguments["--numSequences"]) else "%.5f"%(sum(pos.values())/denominator[p]) ] + [ "%.5f"%(pos.get(r,0)/sum(pos.values())) if sum(pos.values()) > 0 else "0.00" for r in aa_list ] )
	    
		#clean up
		for f in glob.glob("%s.*"%tempFile): 
			os.remove(f)

	print( "Successfully built %d/%d profiles for %s using %d sequences!" % ( success, numProfiles, vgene, len(seqs)-len(germList[vgene]) ) )
	return results

	
def main():

	global masterList, germList

	#load sequences
	masterList = defaultdict( list )
	with open(arguments["<sequences.fa>"], 'r') as handle:
		for sequence in SeqIO.parse(handle, "fasta"):

			#start with a special case where IMGT allele is misnamed
			sequence.description = re.sub("IGHV4-4\*0[78]", "IGHV4-59*11", sequence.description)

			#collapse distal alleles
			sequence.description = re.sub("(V\d-\d+)D", r'\1', sequence.description)

			gene = re.search("(IG[HKL]V\d-[^*]+)", sequence.description) #breaks all the ORs - don't care
			if gene:
				if not arguments["-a"]:
					sequence.seq = sequence.seq.translate() #don't care about anything else in this script
				masterList[ gene.group(1) ].append( sequence )

	#replace list with array
	#have to build manually because otherwise numpy is turning SeqRecords
	#	     into lists of chars (AAs), which causes random.choice
	#	     to throw an error (non 1-D array)
	# (weird footnote: this _doesn't_ happen if 1 or more sequences in the
	#	     don't have an even number of codons)
	# Anyway, only fix I can come up with is to manually place each SeqRecord
	#	     in the array. We have to do it here, afterward, because until
	#	     we've finished loading the sequences, I don't know how many
	#	     there will be of each germline and numpy arrays have to be
	#	     pre-allocated.

	for v in masterList.keys():
		a = numpy.empty( len(masterList[v]), dtype=object )
		for i in range(len(masterList[v])):
			a[i] = masterList[v][i]
		masterList[v] = a


	#load germlines
	germList = defaultdict( list )
	with open(arguments["--germline"], 'r') as handle:
		for sequence in SeqIO.parse(handle, "fasta"):
			#start with a special case where IMGT allele is misnamed
			sequence.id = re.sub("IGHV4-4\*0[78]", "IGHV4-59*11", sequence.id)
			#collapse distal alleles and remove allele designation
			gene = re.sub("(V\d-\d+)D?\*.*", r'\1', sequence.id)
			germList[ gene ].append( sequence )

			if gene not in mask:
				mask[gene] = arguments['--mask']
    

	#now let's start building profiles
	gsspPool = Pool( arguments['-t'] )
	profiles = gsspPool.map( buildGSSP, sorted(masterList.keys()) )
	gsspPool.close()
	gsspPool.join()
	

	#save output
	with open(arguments["--output"], "w") as outHandle:
		output = csv.writer(outHandle, delimiter="\t", dialect='unix', quoting=csv.QUOTE_NONE)
		output.writerow( ["Vgene", "prof#", "pos", "germ", "freq"] + aa_list )
		for blob in profiles:
			for row in blob:
				output.writerow( row )


if __name__ == '__main__':

	aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

	arguments = docopt(__doc__)
	arguments['--numSequences'] = int(arguments['--numSequences'])
	arguments['--profiles'] = int(arguments['--profiles'])
	arguments['-t'] = int(arguments['-t'])

	if arguments['--germline'] == "SONAR/germDB/IgHKLV_cysTruncated.AA.fa":
		arguments['--germline'] = re.sub("SONAR", SCRIPT_FOLDER, arguments['--germline'])

	mask = dict()
	try:
		arguments['--mask'] = int( arguments['--mask'] )
	except ValueError:
		with open( arguments['--mask'], "r" ) as handle:
			reader = csv.reader( handle, delimiter="\t" )
			for row in reader:
				mask[ row[0] ] = int( row[1] )
		arguments['--mask'] = 0 #for all other genes
	

	#log command line
	logCmdLine(sys.argv)
    
	prj_tree = ProjectFolders(os.getcwd())
	os.makedirs("%s/work/mGSSP"%prj_tree.home, exist_ok=True)

	main()

