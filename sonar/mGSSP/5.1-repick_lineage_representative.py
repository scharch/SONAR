#!/usr/bin/env python3

"""
5.1-repick_lineage_representative.py

This script uses the pseudo-lineages identified by 2.4-cluster_into_groups.py.
      It aligns each lineage using the full V(D)J sequence and selects a new
      lineage representative based on the sequences which is closest to the
      lineage consensus.

Usage: 5.1-repick_lineage_representative.py [ -i input.fa -o output.fa -m 2 -t 1 ]

Options:
    -i input.fa    File with lineage-anotated sequences to align.
                      [default: output/sequences/nucleotide/<project>_goodVJ_unique_lineageNotations.fa]
    -o output.fa   Where to save re-picked lineage representatives.
                      [default: output/sequences/nucleotide/<project>_consensusLineageRepresentatives.fa]
    -m 2           Minimum number of sequences in psuedo-lineage to keep. [default: 2]
    -t 1	   Number of threads to use. [default: 1]


Created by Chaim A Schramm on 2016-05-24.
Added to SONAR as part of mGSSP on 2017-02-24.
Modified to use VSearch instead of USearch by CAS on 2018-07-30.
Edited to use Py3 and DocOpt by CAS 2018-08-29.
Multithreaded clustering and alignment by CAS 2018-09-05.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
			 Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
import itertools
from collections import OrderedDict
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline
from multiprocessing import Pool


try:
	from sonar.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/mGSSP")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *


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


def getHandle( item ):
    
	if item in handleList:
		return handleList[item]
	else:
		if len(handleList) >= 250: #being conservative here; increases runtime but decreases possibility of error
			i,h = handleList.popitem(False)
			h.close()
		handleList[item] = open("%s/%s.fa"%(prj_tree.lineage, item), "a")
		return handleList[item]


def closestToConsensus(linIt):

	results = []
	print( "Starting vsearch on a new chunk..." )

	for lineage in linIt:
		#save time on singletons (if they weren't excluded by minSeq)
		if len(lineage['desc']) == 1:
			with open("%s/%s.fa" % (prj_tree.lineage, lineage['name']), "rU") as handle:
				results.append( SeqIO.read(handle,'fasta') )
		else:

			FNULL = open(os.devnull, 'w') #don't clutter up output with tons of vsearch messages

			#cluster and rapid align with vsearch
			subprocess.call([usearch, "-cluster_size", "%s/%s.fa" % (prj_tree.lineage, lineage['name']), 
					 "-id", "0.97", "-sizein", "-sizeout",
					 "-msaout", "%s/%s_msa.fa"%(prj_tree.lineage, lineage['name']), "-clusterout_sort"],
					stdout=FNULL, stderr=subprocess.STDOUT)

			#extract biggest cluster
			with open("%s/%s_msa.fa"%(prj_tree.lineage, lineage['name']), "rU") as allClusters:
				with open("%s/%s_msaBiggest.fa"%(prj_tree.lineage, lineage['name']), "w") as biggestOnly:
					blank = next( allClusters )
					for line in allClusters:
						if "consensus" in line:
							break
						biggestOnly.write(line)
	    
			#open the msa
			with open("%s/%s_msaBiggest.fa" % (prj_tree.lineage, lineage['name']), "rU") as handle:
				aln = AlignIO.read(handle, "fasta")
		
			#add derep size to alignment as weighting
			for rec in aln:
				rec.annotations['weight'] = int( rec.id.split(";")[1].split("=")[1] )

			summary_align = AlignInfo.SummaryInfo(aln)
			pssm = summary_align.pos_specific_score_matrix()

			#score each sequence and save the best one
			scores = dict()
			for record in aln:
				myScore = 0
				for i,l in enumerate(record):
					myScore += pssm[i][l]
				scores[record.id] = myScore
				
			d=sorted(aln, key=lambda rec: scores[rec.id], reverse=True) #reverse -> get max
			d[0].seq = d[0].seq.ungap("-") #remove gaps
			d[0].id = d[0].id.split(";")[0] #remove usearch size annotation
			d[0].id = re.sub("^\*","",d[0].id) #get rid of possible annotation from vsearch
			d[0].description = lineage['desc'][d[0].id] #restore original info

			results.append( d[0] )

	return results


def main():

	#first, open the input file and parse into pseudo-lineages
	lineages = dict()
	count = 0
	for sequence in SeqIO.parse(open(arguments['-i'], "rU"), "fasta"):
		info = re.search(" size=(\d+) lineage_num=(\d+).*lineage_size=(\d+)", sequence.description)
		if info:
			if int(info.group(3)) >= arguments['-m']:
				if info.group(2) not in lineages:
					lineages[info.group(2)] = dict( name=info.group(2), size=int(info.group(3)), desc=dict() )

				lineages[info.group(2)]['desc'][sequence.id] = sequence.description #because otherwise the alignment loses it

				#derep 'singletons' (seq size and lin size are the same) don't need to 
				#    be renamed for usearch and don't need saved handles
				if info.group(1) == info.group(3):
					with open("%s/%s.fa"%(prj_tree.lineage ,info.group(2)), "w") as handle:
						SeqIO.write([sequence], handle, 'fasta')
				else:
					#add size notation in a way usearch can understand
					sequence.id += ";size=%s" % info.group(1) #it's a string because that's what re.search returns
					SeqIO.write([sequence], getHandle(info.group(2)), 'fasta')

		count += 1
		if count % 100000 == 0: print( "Processed %d sequences in %d lineages so far..." % (count, len(lineages)) )

	#cleanup
	for h in handleList.values(): h.close()
	handleList.clear()

	
	#go through each lineage and do the alignment
	reps = list()
	
	if arguments['-t'] > 1:
		pool = Pool(arguments['-t'])
		blob = pool.map( closestToConsensus, iterator_slice(lineages.values(), 1000) )
		pool.close()
		pool.join()

		for result in blob:
			reps += result

	else:
		#don't thread
		reps = closestToConsensus( lineages.values() )

	    
	#write output
	with open( arguments['-o'], "w" ) as handle:
		SeqIO.write( reps, handle, "fasta" )




if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['-m'] = int( arguments['-m'] )
	arguments['-t'] = int( arguments['-t'] )

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	arguments['-i'] = re.sub("<project>", prj_name, arguments['-i'])
	arguments['-o'] = re.sub("<project>", prj_name, arguments['-o'])

	if not os.path.isfile(arguments['-i']):
		sys.exit( "Can't find input file %s" % arguments['-i'] )

	       
	#log command line
	logCmdLine(sys.argv)	

	
	#since we're having to open/close/reopen filehandles to manage system limits
	#  clear out any old files so we're not double-counting data from multiple runs
	print("Cleaning up old files...", end="", flush=True)
	for f in glob.glob("%s/*fa"%prj_tree.lineage):
		os.remove(f)
	print("Done")
	#make dir for msas but ignore error if it already exists
	os.makedirs("%s/msa"%prj_tree.lineage, exist_ok=True)
	
	handleList = OrderedDict()

	
	main()

