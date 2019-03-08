#!/usr/bin/env python3

"""
1.4-cluster_sequences.py

This script is a mostly a simple wrapper for vsearch that performs two steps 
     of clustering to group sequences that may differ only due to sequencing 
     error. The first step finds exact duplicates, and potentially removes 
     reads with coverage lower than a specified threshold. The second step 
     further clusters the dereplicated sequences at a lower sequence 
     identity, preferentially using reads with higher coverage as centroids.
Output stats will also be added to the AIRR-formatted rearrangements file.

Usage: 1.4-cluster_sequences.py [ -f input.fa --min1 1 --min2 3 --id .99 --maxgaps 0 ]

Options:
    -f input.fa   Input sequence file in fasta format. 
                     [default: output/sequences/nucleotide/<project>_goodVJ.fa]
    --min1 1      Minimum identical replicates needed to keep a read in the first
                     step of clustering. [default: 1]
    --min2 3      Minimum number of total reads (ie counting exact replicates
                     from step 1 separately) in a cluster needed to keep a
                     sequence in the second step of clustering. [default: 3]
    --id .99      Percent sequence identity used for the second step of
                     clustering. [default: 0.99]
    --maxgaps 0   vsearch parameter specifying how many (non-terminal) gaps are
                     allowed in an alignment for two sequences to cluster together.
                     Should not be changed, in most cases. [default: 0]

Created by Zizhang Sheng.
Edited to switch to VSearch by Chaim A Schramm 2018-07-30.
Ported to Python to handle AIRR-format data by CAS 2018-10-16.
Added self-reference for centroid sequences and deprecated use of 'unique' status
    by CA Schramm 2019-03-07.
Added manual maxgaps option by CA Schramm 2019-03-08.

Copyright (c) 2011-2019 Columbia University and Vaccine Research Center, National 
                         Institutes of Health, USA. All rights reserved.

"""


import sys
from docopt import docopt
from Bio import SeqIO
import airr

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *


	
def reformatInput(seqFile):
	with open(seqFile, "r") as seqs:
		for s in SeqIO.parse(seqs, "fasta"):
			mySize = re.search( "duplicate_count=(\d+)", s.description )
			if mySize:
				s.id += ";size=%s" % mySize.group(1)
			yield s



def getUniques(seqFile, size_dict):
	with open(seqFile, "r") as seqs:
		for s in SeqIO.parse(seqs, "fasta"):
			if s.id in size_dict:
				s.description += " cluster_count=%d" % size_dict[s.id]
				yield s


				
def main():

	#start by making possible "duplicate_count" info available to vsearch
	with open("temp.fa", "w") as handle:
		SeqIO.write( reformatInput(arguments['-f']), handle, "fasta" )
	
	#first step on higher identity
	subprocess.call( [ vsearch, "-derep_fulllength", "temp.fa",
			   "-output", "temp_dedup.fa",
			   "-uc", "temp.uc",
			   "-sizein", "-sizeout",
			   "-minuniquesize", arguments['--min1'] ] )

	#process the uc file
	centroid = dict()
	with open("temp.uc" ,"r") as handle:
		uc = csv.reader(handle, delimiter="\t")
		for row in uc:
			if row[0] == "H":
				centroid[ re.sub(";size=\d+","",row[8]) ] = re.sub(";size=\d+","",row[9])

	#second clustering step
	subprocess.call( [ vsearch, "-cluster_size", "temp_dedup.fa",
			   "-sizein", "-sizeout",
			   "-maxgaps", arguments['--maxgaps'],
			   "-id", arguments['--id'],
			   "-uc", "%s.cluster"%os.path.splitext(arguments['-f'])[0] ] )

	#process the uc file
	size = dict()
	with open("%s.cluster"%os.path.splitext(arguments['-f'])[0], "r") as handle:
		uc = csv.reader(handle, delimiter="\t")
		for row in uc:
			if row[0] == "H":
				centroid[ re.sub(";size=\d+","",row[8]) ] = re.sub(";size=\d+","",row[9])
			elif row[0] == "C":
				#have the centroids point to themselves for more uniform dowsntream processing
				centroid[ re.sub(";size=\d+","",row[8]) ] = re.sub(";size=\d+","",row[8])
				#but only save them if they meet the threshold
				if int(row[2]) >= arguments['--min2']:
					size[ re.sub(";size=\d+","",row[8]) ] = int(row[2])

	#clean up
	os.remove("temp.fa")
	os.remove("temp_dedup.fa")
	os.remove("temp.uc")

	#do sequence outputs
	with open("%s_unique.fa"%os.path.splitext(arguments['-f'])[0], "w") as handle:
		SeqIO.write( getUniques(arguments['-f'], size), handle, 'fasta' )
		  
	#retrieve unique CDR3s (and do AA seqs as appropriate)
	if "goodVJ" in arguments['-f']:
		cdr3_file = re.sub("goodVJ","goodCDR3", arguments['-f'])
		if os.path.isfile(cdr3_file):
			with open("%s_unique.fa"%os.path.splitext(cdr3_file)[0], "w") as handle:
				SeqIO.write( getUniques(cdr3_file, size), handle, 'fasta' )
		else:
			print( "Can't find %s to extract unique sequences..."%cdr3_file, file=sys.stderr )
			
		if "nucleotide" in cdr3_file:
			cdr3_file = re.sub("nucleotide","amino_acid", cdr3_file)
			if os.path.isfile(cdr3_file):
				with open("%s_unique.fa"%os.path.splitext(cdr3_file)[0], "w") as handle:
					SeqIO.write( getUniques(cdr3_file, size), handle, 'fasta' )
			else:
				print( "Can't find %s to extract unique sequences..."%cdr3_file, file=sys.stderr )
	if "nucleotide" in arguments['-f']:
		aa_file = re.sub("nucleotide","amino_acid", arguments['-f'])
		if os.path.isfile(aa_file):
			with open("%s_unique.fa"%os.path.splitext(aa_file)[0], "w") as handle:
				SeqIO.write( getUniques(aa_file, size), handle, 'fasta' )
		else:
			print( "Can't find %s to extract unique sequences..."%aa_file, file=sys.stderr )

	#now do AIRR output
	if "output/sequences/nucleotide" in arguments['-f']:
		if os.path.isfile("%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name)):
			clustered = airr.derive_rearrangement( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name),
							       fields=['centroid', 'cluster_count'] )
			for r in airr.read_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) ):
				#clear old annotations in case we ran 1.4 previously
				r['centroid'] = ""
				r['cluster_count'] = ""

				#now add back current annotations
				#two rounds of clustering means we start by looking for the centroid of the centroid,
				#    falling back to the first level centroid (second clustering step only) if appropriate
				r['centroid'] = centroid.get( centroid.get(r['sequence_id'],""), centroid.get(r['sequence_id'],"") )

				#add cluster size information for final centroids. I am doing away with changing the 'status' of
				#    the centroids to 'unique' because I've started using this script in a lot of cases where it
				#    doesn't make sense to treat 'unique' as a subset of 'good', and I therefore need to preserve
				#    the original status designation. To find centroids, look for a non-null 'cluster_count' field
				if r['sequence_id'] in size:
					r['cluster_count'] = size[ r['sequence_id'] ]

				clustered.write(r)
			clustered.close()
			os.rename( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) )
		else:
			print( "Can't find the rearrangements file, not saving data in AIRR format", file=sys.stderr )



if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--min2'] = int( arguments['--min2'] )

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	arguments['-f'] = re.sub("<project>", prj_name, arguments['-f'])

	if not os.path.isfile(arguments['-f']):
		sys.exit( "Can't find input file %s" % arguments['-f'] )

	       
	#log command line
	logCmdLine(sys.argv)

	
	main()

