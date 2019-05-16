#!/usr/bin/env python3

"""
2.3-intradonor_analysis.py

This script uses an iterative phylogenetic analysis to find sequences related
      to a set of known antibodies. Preprocessed sequences (including optional
      filtering by assigned germline) are split into groups and used together
      with known sequences to build neighbor-joining trees rooted on the
      germline V gene of the known antibodies. Sequences in the minimum
      sub-tree spanning all of the known sequences are passed forward into the
      next iteration. The algorithm is considered to have converged when 95%
      of the input sequences in a round are in the minimum sub-tree.

This script has an option to run the analysis as a cluster job. If run locally,
      can be multithreaded with use of the -threads parameter.

This algorithm is generally intended to find somatically related antibodies
      from a single lineage within a single donor. However, in the special
      case of VRC01 class antibodies, we have shown that exogenous VRC01 class
      heavy chains can be used for "cross-donor" analysis to identify a 
      lineage of VRC01 class antibodies within a new donor. See, eg:
      * Wu, Zhou, Zhu et al. Science (2011) doi: 10.1126/science.1207532
      * Zhu et al PNAS (2013) doi: 10.1073/pnas.130626211
      * Wu, Zhang, Schramm, Joyce, Do Kwon, Zhou, Sheng, et al. Cell (2015)
                                               doi: 10.1016/j.cell.2015.03.004
      It is unknown at this time whether this method would work for other
      classes of antibodies.

Usage: 2.3-intradonor_analysis.py --n native.fa --v germline_V [--locus H | --lib path/to/library.fa] [--in custom/input.fa | [ -a --noFilter ]] [options]

Options:
    --n	<native.fa>       Fasta file containing the known sequences.
    --v	<IGHV1-2*02>      Assigned germline V gene of known antibodes, for use in 
                             rooting the trees.
    --locus H             H: use V heavy / K: use V kappa / L: use V lambda
                             PLEASE NOTE: This algorithm is not generally recommended
                             for analyzing light chain sequences.
                             [default: H]
    --lib <germline.fa>   Optional custom germline library (eg Rhesus or Mouse).
                             Mutually exclusive with --locus.
    --in <seqs.fa>        Optional custom set of sequences to be analayzed.
                             Default is output/sequences/<ROOT>_goodVJ_unique.fa (or 
                             output/sequences/<ROOT>_allV.fa with the -a flag).
    -a		          Use all NGS sequences with an assigned V, even those with
                             without an assigned J gene or other possible errors.
                             Not available in combination with --in.
                             Instead of using -a, recommended usage is to run
                             `1.4-dereplicate_sequences.pl -f output/sequences/<ROOT>_allV.fa -s 2`
                             and then calling this script with `--in output/sequences/ROOT_allV_unique.fa`
                             [default: False]
    --noFilter            Do NOT filter NGS sequences for correct germline V gene
                             assignment. Not available in combination with --in
                             (which assumes you've done your own filtering, if you
                             want it). [default: False]
    --cluster             Submit tree-building jobs to the cluster. [default: False]
    --maxIters <15>       Optional maximum number of rounds to conduct before giving up.
                             [default: 15]
    --npf <250>           Optional number of sequences to include in each split file.
                             Larger number results in slower runtime due to constructing
		             the MSA but fewer iterations to convergence. Default = 250
		             when running locally and 1,000 on a cluster.
    --threads <1>         Number of threads to use when running locally. [default: 1]
    -f                    Force a restart of the analysis, even if there are files from
                             a previous run in the working directory. [default: False]

Created by Zhenhai Zhang on 2011-07-12.
Modified to compress into a single script and many updates by 
        Chaim A Schramm 2014-01-09.
Edited and commented for publication by Chaim A Schramm on 2015-04-14.
Edited to use Py3 and DocOpt by CAS 2018-08-22.

Copyright (c) 2011-2018 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import time, sys
from docopt import docopt
from multiprocessing import Pool
from functools import partial
from io import StringIO
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline

try:
	from SONAR.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/lineage")
	sys.path.append(find_SONAR[0])
	from SONAR.lineage import *



def muscleProcess (threadID, filebase, outbase, treebase):

	fasta	 = filebase % threadID
	output	 = outbase  % threadID
	treeFile = treebase % threadID

	print( "Building NJ tree from %s" % fasta )

	run_muscle = MuscleCommandline( cmd=muscle, input=fasta, out=output )
	run_muscle.tree1      = treeFile
	run_muscle.cluster1   = "neighborjoining"
	run_muscle.maxiters   = 1
	thisVarHidesTheOutput = run_muscle()




def main():

	global converged, num_nats, germ_seq
	currentIter = 0
	log = open( "%s/intradonor.log" % prj_tree.logs, "a+" )

	# master loop
	while not converged:
		
		# parse tree files and get all reads clustered with native antibodies
		tree_files, good, total, retained_reads = sorted(glob.glob("%s/NJ*.tree" %prj_tree.lineage)), 0.0, 0.0, []
								 
		if arguments['-f']:
			tree_files = [] #no need to process the files since we will be starting over
			arguments['-f'] = False #turn it off so we don't get stuck in an infinite loop of restarts

		#this gets skipped in the first round
		for idx, tf in enumerate(tree_files):
			
			#read in the tree as a string
			tree_string = open(tf, "r").read().strip()
			
			#get rid of possible negative branch lengths and then parse
			tree_string, num = re.subn(":-\d+\.\d+",":0", tree_string)
			tree_string, num = re.subn("\n","", tree_string)

			#need to supply comments_are_confidence=True to prevent Phylo from interpreting numeric sequence IDs as bootstrap support
			tree = Phylo.read(StringIO(tree_string), "newick", comments_are_confidence=True)

			
			#outgroup root the tree on the germline
			try:
				germ_node = next( tree.find_clades(arguments['--v']) )
			except StopIteration:
				sys.exit( "Can't find germline gene %s in file %s" % (arguments['--v'], tf) )
			tree.root_with_outgroup(germ_node)


			# start by finding one of the natives
			try:
				nat1 = next( tree.find_clades(natives_list[0]) )
			except StopIteration:
				sys.exit( "Can't find native antibody %s in file %s" % (natives_list[0], tf) )

				
			#walk up the tree from this starting point until we find the minimum sub-tree with all natives
			path_to_nat1 = tree.get_path(nat1)
			level = -2
			subtree = path_to_nat1[level]
			all_leaves = [node.name for node in subtree.get_terminals()]
			try:
				while not all( [ nat in all_leaves for nat in natives_list ] ):
					level -= 1
					subtree = path_to_nat1[level]
					all_leaves = [node.name for node in subtree.get_terminals()]
			except IndexError:
				sys.exit( "Can't find a subtree with all native sequences in file %s" % tf )

			#check if germline sequence is in the subtree to keep counts correct
			if arguments['--v'] in all_leaves:
				all_leaves.remove(arguments['--v'])
				
			# save sequences in subtree and print progress message
			retained_reads += all_leaves
			good += len(all_leaves) - num_nats
			total += len(tree.get_terminals()) - num_nats - 1 #also don't count germline
			if not arguments['--cluster']:
				print( "Found %d reads in subtree #%d. Total saved so far: %d / %d" % ( len(all_leaves)-num_nats, idx+1, good, total-1) )


		# processed all trees from last round, now do a sanity check
		if good == 0 and len(tree_files) > 0:
			log.write( "%s - Round %d: NO positive sequences found --stopped!\n" % (time.strftime("%H:%M:%S"), currentIter) )
			log.close()
			sys.exit( "NO positive sequences found --stopped!" )

		
		# Are we starting a new run? If not, limit memory usage by only reading in the retained reads
		read_dict = dict()
		if len(tree_files) == 0:
			log.write( "%s - Starting a new analysis from scratch...\n" % time.strftime("%H:%M:%S") )
			print( "%s - Starting a new analysis from scratch..." % time.strftime("%H:%M:%S") )
			if not arguments['--noFilter']:
				#load by gene but ignore allele
				read_dict = load_fastas_with_Vgene( arguments['--in'], arguments['--v'].split("*")[0] )
			else:
				read_dict = load_fastas( arguments['--in'] )
		else:
			read_dict = load_seqs_in_dict( arguments['--in'], set(retained_reads) )

		#error checking
		if len(read_dict) == 0:
			log.write( "%s - Error: failed to load any sequences from %s, stopped\n" % (time.strftime("%H:%M:%S"), arguments['--in']) )
			log.close()
			sys.exit( "Error: failed to load any sequences from %s, stopped" % arguments['--in'] )

		#randomize the order
		shuffled_reads = list(read_dict.values())
		random.shuffle(shuffled_reads)
			
		# Check for convergence before starting a new round
		if total == 0 or good/total < 0.95: #total == 0 would be round 1, so don't want to quit early
			
			if currentIter >= arguments['--maxIters']:
				log.write( "%s - Maximum number of iterations reached without convergence. Current round: %d reads, %5.2f%% of input\n" % (time.strftime("%H:%M:%S"), good, 100*good/total) )
				log.close()
				sys.exit( "Maximum number of iterations reached without convergence. Current round: %d reads, %5.2f%% of input" % (good, 100*good/total) )
			else:
				if len(tree_files) > 0: 
					#it's a silly message to print the first time through
					log.write( "%s - Finished processing round %d: %d reads, %5.2f%% of input\n" % (time.strftime("%H:%M:%S"), currentIter, good, 100*good/total) )
					print( "%s - Finished processing round %d: %d reads, %5.2f%% of input" % (time.strftime("%H:%M:%S"), currentIter, good, 100*good/total) )
				currentIter +=1


			#do some cleanup
			oldFiles = glob.glob("%s/*" % prj_tree.last)
			for old in oldFiles:
				os.remove(old)
			lastRound = glob.glob("%s/NJ*" % prj_tree.lineage)
			for infile in lastRound:
				os.rename( infile, "%s/%s" % (prj_tree.last, os.path.basename(infile)) )
					
			#open initial output
			f_ind = 1
			outFasta = open("%s/NJ%05d.fa" % (prj_tree.lineage, f_ind), "w")
			currentSize = 0
			
			#process reads
			for r in shuffled_reads:

				outFasta.write(">%s\n%s\n" % (r.id, r.seq))
				currentSize += 1
					
				#reached max size for this split?
				if currentSize == arguments['--npf']:

					#add native and germline sequences before closing file
					outFasta.write( ">%s\n%s\n" % (germ_seq.id, germ_seq.seq) )
					for n in natives.values():
						outFasta.write( ">%s\n%s\n" % (n.id, n.seq) )
					outFasta.close()
					f_ind += 1
					outFasta = open("%s/NJ%05d.fa" % (prj_tree.lineage, f_ind), "w")
					currentSize = 0

			#Any leftovers?
			if currentSize > 0:
				outFasta.write( ">%s\n%s\n" % (germ_seq.id, germ_seq.seq) )
				for n in natives.values():
					outFasta.write( ">%s\n%s\n" % (n.id, n.seq) )
				outFasta.close()
			else:
				outFasta.close()
				f_ind -= 1 #don't submit an empty file to the cluster

				
			# At this point, if we are running on the cluster, submit this round and quit loop
			if arguments['--cluster']:
				# write pbs files and auto submit shell script
				command = "NUM=`printf \"%%05d\" $SGE_TASK_ID`\n%s -in %s/NJ$NUM.fa -out %s/NJ$NUM.aln -cluster1 neighborjoining -maxiters 1 -tree1 %s/NJ$NUM.tree" % \
				    (muscle, prj_tree.lineage, prj_tree.lineage, prj_tree.lineage)
				pbs = open("%s/intradonor.sh" % prj_tree.lineage, 'w')
				pbs.write( PBS_STRING%("%s-intradonor"%prj_name, "500M", "1:00:00", "%s > %s/NJ$NUM.out 2> %s/NJ$NUM.err"%(command, prj_tree.lineage, prj_tree.lineage)) )
				pbs.close()

				# write a second PBS job to restart this script once muscle has finished
				next_round = open("%s/nextround.sh" % prj_tree.lineage, 'w')
				next_round.write( "\
#!/bin/bash\n\
#$ -hold_jid %s-intradonor\t\t# run once previous round has completed\n\
#$ -N restart-intradonor\t\t# job name\n\
#$ -l mem=1G,time=1:00:00\t\t# resource requests\n\
#$ -cwd\t\t\t\t\t# use current directory\n\
#$ -o %s/restart-intradonor.out\t#output\n\
#$ -e %s/restart-intradonor.err\t#error\n\
%s/lineage/2.3-intradonor_analysis.py --n %s --v %s --lib %s --in %s --maxIters %d --cluster\n" % 
						  (prj_name, prj_tree.lineage, prj_tree.lineage, SCRIPT_FOLDER, arguments['--n'], arguments['--v'], arguments['--library'], arguments['--in'], arguments['--maxIters']-1) )
				next_round.close()

				os.system( "%s -t 1-%d %s/intradonor.sh" % (qsub, f_ind, prj_tree.lineage) )
				os.system( "%s %s/nextround.sh"	 % (qsub, prj_tree.lineage) )
				log.write( "%s - Submitted current round to cluster\n" % time.strftime("%H:%M:%S") )
				log.close()
				break #exits "while not converged" loop

			else:
				#run locally
				partial_muscle = partial( muscleProcess, filebase="%s/NJ%%05d.fa"%prj_tree.lineage, 
							  outbase="%s/NJ%%05d.aln"%prj_tree.lineage, treebase="%s/NJ%%05d.tree"%prj_tree.lineage )
				muscle_pool = Pool(arguments['--threads'])
				muscle_pool.map(partial_muscle, range(1,f_ind+1))
				muscle_pool.close()
				muscle_pool.join()
				

		
		# If we have converged, write final output
		else:
			out_nt	 = open("%s/%s_intradonor_positives.fa"	 % (prj_tree.nt, prj_name),	"w")
			out_aa	 = open("%s/%s_intradonor_positives.fa"	 % (prj_tree.aa, prj_name),	"w")
			SeqIO.write(read_dict.values(), out_nt, "fasta")
			SeqIO.write( [SeqRecord(r.seq.translate(),id=r.id, description=r.description) for r in read_dict.values()], out_aa, "fasta" )
			out_nt.close()
			out_aa.close()

			out_list = open("%s/%s_intradonor_positives.txt" % (prj_tree.tables, prj_name), "w")
			for r in read_dict.values():
				out_list.write( "%s\n" % r.id )
			out_list.close()
			
			log.write( "%s - After round #%d: tree has converged with %d reads!\n" % (time.strftime("%H:%M:%S"), currentIter, good) )
			log.close()
			print( "After round #%d: tree has converged with %d reads!" % (currentIter,good) )
			converged = True

	

if __name__ == '__main__':

	arguments = docopt(__doc__)


	#load native sequences
	if not os.path.isfile(arguments['--n']):
		print("Error: cannot find file with native sequences")
		sys.exit(1)
	natives	      =	 load_fastas(arguments['--n'])
	natives_list  =	 list( natives.keys() )
	num_nats      =	 len(natives_list)

	
	converged = False #keep track of convergence when running locally
	
	
	arguments['--maxIters'] = int( arguments['--maxIters'] )
	arguments['--threads']	= int( arguments['--threads'] )
	
	if arguments['--npf'] is not None:
		arguments['--npf'] = int( arguments['--npf'] )

	if arguments['--cluster']:
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")
		if arguments['--npf'] is None:
			arguments['--npf'] = 1000
	else:
		if arguments['--npf'] is None:
			arguments['--npf'] = 250
	      

	#get input
	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)
	if arguments['--in'] is not None:
		arguments['--noFilter'] = False
		if not os.path.isfile(arguments['--in']):
			print( "Can't find input file!" )
			sys.exit(1)
	else:
		if arguments['-a']:
			arguments['--in'] = "%s/%s_allV.fa" % (prj_tree.nt, prj_name)
		else:
			arguments['--in'] = "%s/%s_goodVJ_unique.fa" % (prj_tree.nt, prj_name)

	#get germline
	if arguments['--lib'] is not None:
		if not os.path.isfile(arguments['--lib']):
			print( "Can't find custom V gene library file!" )
			sys.exit(1)		   
	else:
		if arguments['--locus'] in dict_vgerm_db.keys():
			arguments['--lib'] = dict_vgerm_db[arguments['--locus']]
		else:
			print("Error: valid options for --locus are H, K, L, KL, and HKL only")
			sys.exit(1)

	germ_dict = load_fastas(arguments['--lib'])
	try:
		germ_seq = germ_dict[ arguments['--v'] ]
	except:
		print( "Specified germline gene (%s) is not present in the %s library!\n" % (arguments['--v'], arguments['--lib']) )
		sys.exit(1)

	#create necessary dirs to avoid errors
	os.makedirs(prj_tree.aa, exist_ok=True)
	os.makedirs(prj_tree.nt, exist_ok=True)
	os.makedirs(prj_tree.lineage, exist_ok=True)
	os.makedirs(prj_tree.last, exist_ok=True)
	os.makedirs(prj_tree.logs, exist_ok=True)
	os.makedirs(prj_tree.tables, exist_ok=True)

	
	#log command line
	logCmdLine(sys.argv)


	main()

