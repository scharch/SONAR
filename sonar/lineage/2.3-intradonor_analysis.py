#!/usr/bin/env python

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

Usage: 2.3-intradonor_analysis.py -n native.fa -v germline_V
                                 [-locus <H|K|L|C> -lib path/to/library.fa
				  -i custom/input.fa -maxIters 15
				  -cluster -npf 250 -threads 1
				  -nofilter -a -h -f]

    Invoke with -h or --help to print this documentation.

    Required parameters:
    n		Fasta file containing the known sequences.
    v		Assigned germline V gene of known antibodes, for use in 
                   rooting the trees.

    Optional parameters:	   
    locus	H (default): use V heavy / K: use V kappa / L: use V lambda /
                   C: use a custom germline V library (specify -lib).
		   PLEASE NOTE: This algorithm is not generally recommended
		   for analyzing light chain sequences.
    lib		Optional custom germline library (eg Rhesus or Mouse).
                   Ignored unless "-locus C" is used.
    i		Optional custom set of sequences to be analayzed. "-nofilter"
                   and "-a" flags will be ignored if this option is used.
		   Default is output/sequences/ROOT_goodVJ_unique.fa (or 
		   output/sequences/ROOT_allV.fa with the -a flag below)
    maxIters	Optional maximum number of rounds to conduct before giving up.
                   Default = 15.
    npf         Optional number of sequences to include in each split file.
                   Larger number results in slower runtime due to constructing
		   the MSA but fewer iterations to convergence. Default = 250
		   when running locally and 1,000 on a cluster.
    threads     Number of threads to use when running locally. Default = 1.

    Optional flags:
    cluster	Submit tree-building jobs to the cluster.
    nofilter	Do NOT filter NGS sequences for correct germline V gene
                   assignment. Default = OFF (DO filter).
    a		Use all NGS sequences with an assigned V, even those with
                   out-of frame junctions and/or stop codons or without a 
		   successfuly assigned J gene. Default = OFF (use in-frame ORF
		   sequences only). Instead of using -a, recommended usage is 
		   to manually dereplicate the output/sequences/ROOT_allV.fa
		   file using 1.4-dereplicate_sequences.pl with the -f flag
		   and then passing that output to the -i parameter of this 
		   script.
    f		Force a restart of the analysis, even if there are files from
                   a previous run in the working directory.

Created by Zhenhai Zhang on 2011-07-12.
Modified to compress into a single script and many updates by 
        Chaim A Schramm 2014-01-09.
Edited and commented for publication by Chaim A Schramm on 2015-04-14.

Copyright (c) 2011-2016 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import time, sys, threading
from cStringIO import StringIO
from sonar.lineage import *
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline

try:
	from sonar.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/lineage")
	sys.path.append(find_SONAR[0])
	from sonar.lineage import *



class muscleThread (threading.Thread):
	def __init__(self, threadID, fasta, output, treeFile):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.fasta    = fasta
		self.output   = output
		self.tree     = treeFile	
	def run(self):
		run_muscle = MuscleCommandline( input=self.fasta, out=self.output )
		run_muscle.tree1      = self.tree
		run_muscle.cluster1   = "neighborjoining"
		run_muscle.maxiters   = 1
		thisVarHidesTheOutput = run_muscle()




def main():

	global npf, converged, maxIters, cluster, force, num_nats, correct_V_only, germlineV, germ_seq, inFile, natFile, locus, library
	currentIter = 0
	log = open( "%s/intradonor.log" % prj_tree.logs, "a+" )

	# master loop
	while not converged:
		
		# parse tree files and get all reads clustered with native antibodies
		# start total at 1 so good/total < .95 when we start a new analysis
		tree_files, good, total, retained_reads = sorted(glob.glob("%s/NJ*.tree" %prj_tree.lineage)), 0.0, 1.0, []
								 
		if force:
			tree_files = [] #no need to process the files since we will be starting over
			force = False #turn it off so we don't get stuck in an infinite loop of restarts

	        #this gets skipped in the first round
		for idx, tf in enumerate(tree_files):
			
		        #read in the tree as a string
			tree_string = open(tf, "rU").read().strip()
			
			#get rid of possible negative branch lengths and then parse
			tree_string, num = re.subn(":-\d+\.\d+",":0", tree_string)
			tree_string, num = re.subn("\n","", tree_string)

			#need to supply comments_are_confidence=True to prevent Phylo from interpreting numeric sequence IDs as bootstrap support
			tree = Phylo.read(StringIO(tree_string), "newick", comments_are_confidence=True)

			
			#outgroup root the tree on the germline
			try:
				germ_node = tree.find_clades(germlineV).next()
			except StopIteration:
				sys.exit( "Can't find germline gene %s in file %s" % (germlineV, tf) )
			tree.root_with_outgroup(germ_node)


			# start by finding one of the natives
			try:
				nat1 = tree.find_clades(natives_list[0]).next()
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
		
				
			# save sequences in subtree and print progress message
			retained_reads += all_leaves
			good += len(all_leaves) - num_nats
			total += len(tree.get_terminals()) - num_nats - 1 #also don't count germline
			if idx % 25 == 0 and not cluster:
				print "Found %d reads in subtree #%d. Total saved so far: %d / %d" % ( len(all_leaves)-num_nats, idx+1, good, total-1)


		# processed all trees from last round, now do a sanity check
		if good == 0 and len(tree_files) > 0:
			log.write( "%s - Round %d: NO positive sequences found --stopped!\n" % (time.strftime("%H:%M:%S"), currentIter) )
			log.close()
			sys.exit( "NO positive sequences found --stopped!" )

		
		# Are we starting a new run? If not, limit memory usage by only reading in the retained reads
		read_dict = dict()
		if len(tree_files) == 0:
			log.write( "%s - Starting a new analysis from scratch...\n" % time.strftime("%H:%M:%S") )
			print "%s - Starting a new analysis from scratch..." % time.strftime("%H:%M:%S")
			if correct_V_only:
				#load by gene but ignore allele
				read_dict = load_fastas_with_Vgene( inFile, germlineV.split("*")[0] )
			else:
				read_dict = load_fastas( inFile )
		else:
			read_dict = load_seqs_in_dict( inFile, set(retained_reads) )

		#error checking
		if len(read_dict) == 0:
			log.write( "%s - Error: failed to load any sequences from %s, stopped\n" % (time.strftime("%H:%M:%S"), inFile) )
			log.close()
			sys.exit( "Error: failed to load any sequences from %s, stopped" % inFile )

		#randomize the order
		shuffled_reads = read_dict.values()
		random.shuffle(shuffled_reads)
			
		# Check for convergence before starting a new round
		if good/total < 0.95:
			
			if currentIter >= maxIters:
				log.write( "%s - Maximum number of iterations reached without convergence. Current round: %d reads, %5.2f%% of input\n" % (time.strftime("%H:%M:%S"), good, 100*good/total) )
				log.close()
				sys.exit( "Maximum number of iterations reached without convergence. Current round: %d reads, %5.2f%% of input" % (good, 100*good/total) )
			else:
				if len(tree_files) > 0: 
					#it's a silly message to print the first time through
					log.write( "%s - Finished processing round %d: %d reads, %5.2f%% of input\n" % (time.strftime("%H:%M:%S"), currentIter, good, 100*good/total) )
					print "%s - Finished processing round %d: %d reads, %5.2f%% of input" % (time.strftime("%H:%M:%S"), currentIter, good, 100*good/total)
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
				if currentSize == npf:

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
			if cluster:
				# write pbs files and auto submit shell script
				command = "NUM=`printf \"%%05d\" $SGE_TASK_ID`\n%s -in %s/NJ$NUM.fa -out %s/NJ$NUM.aln -cluster1 neighborjoining -maxiters 1 -tree1 %s/NJ$NUM.tree" % \
				    (cluster_muscle, prj_tree.lineage, prj_tree.lineage, prj_tree.lineage)
				pbs = open("%s/intradonor.sh" % prj_tree.lineage, 'w')
				pbs.write( PBS_STRING%(f_ind, "%s-intradonor"%prj_name, "500M", "1:00:00", "%s > %s/NJ$NUM.out 2> %s/NJ$NUM.err"%(command, prj_tree.lineage, prj_tree.lineage)) )
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
%s/lineage/2.3-intradonor_analysis.py -n %s -v %s -locus %s -lib %s -i %s -maxIters %d -cluster\n" % 
						  (prj_name, prj_tree.lineage, prj_tree.lineage, SCRIPT_FOLDER, natFile, germlineV, locus, library, inFile, maxIters-1) )
				next_round.close()

				os.system("qsub %s/intradonor.sh" % prj_tree.lineage)
				os.system("qsub %s/nextround.sh"  % prj_tree.lineage)
				log.write("%s - Submitted current round to cluster\n" % time.strftime("%H:%M:%S"))
				log.close()
				break #exits "while not converged" loop

			else:
				#run locally
				currentFile = 1
				allThreads = []
				while currentFile <= f_ind:
					if threading.activeCount() <= numThreads:
						muscle = muscleThread( currentFile, "%s/NJ%05d.fa" % (prj_tree.lineage, currentFile),
								       "%s/NJ%05d.aln" % (prj_tree.lineage, currentFile),
								       "%s/NJ%05d.tree" % (prj_tree.lineage, currentFile) )
						print "Building NJ tree from %s/NJ%05d.fa" % (prj_tree.lineage, currentFile)
						muscle.start()
						allThreads.append(muscle)
						currentFile += 1
					else:
						#queue is full, so take a break
						time.sleep(60)
				for t in allThreads:
					t.join()


		
		# If we have converged, write final output
		else:
			out_nt   = open("%s/%s_intradonor_positives.fa"  % (prj_tree.nt, prj_name),     "w")
			out_aa   = open("%s/%s_intradonor_positives.fa"  % (prj_tree.aa, prj_name),     "w")
			SeqIO.write(read_dict.values(), out_nt, "fasta")
			SeqIO.write( [SeqRecord(r.seq.translate(),id=r.id, description=r.description) for r in read_dict.values()], out_aa, "fasta" )
			out_nt.close()
			out_aa.close()

			out_list = open("%s/%s_intradonor_positives.txt" % (prj_tree.tables, prj_name), "w")
			for r in read_dict.values():
				out_list.write( "%s\n" % r.id )
			out_list.close()
			
			log.write( "%s - Tree has converged with %d reads!\n" % (time.strftime("%H:%M:%S"), good) )
			log.close()
			print "Tree has converged with %d reads!" % good
			converged = True

	

if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)


	#check forcing parameter
	force = False
	flag = [x for x in ["f", "-f", "--f", "force", "-force", "--force"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		force = True


	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	= fullpath2last_folder(prj_tree.home)


	#check filtering parameter
	correct_V_only = True
	flag = [x for x in ["nof", "-nof", "--nof", "nofilter", "-nofilter", "--nofilter"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		correct_V_only = False

	#check sequence quality parameter
	selectedFile = "%s/%s_goodVJ_unique.fa" % (prj_tree.nt, prj_name)
	flag = [x for x in ["a", "-a", "--a", "all", "-all", "--all"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		selectedFile = "%s/%s_allV.fa" % (prj_tree.nt, prj_name)

	#check cluster parameter
	cluster = False
	default_npf = 250
	flag = [x for x in ["c", "-c", "--c", "cluster", "-cluster", "--cluster"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		cluster = True
		default_npf = 1000

	#get the parameters from the command line
	dict_args = processParas(sys.argv, n="natFile", v="germlineV", locus="locus", lib="library", i="inFile", maxIters="maxIters", npf="npf", threads="numThreads")
	defaults = dict(locus="H", library="", inFile=selectedFile, maxIters=15, npf=default_npf, numThreads=1)
	natFile, germlineV, locus, library, inFile, maxIters, npf, numThreads = getParasWithDefaults(dict_args, defaults, "natFile", "germlineV", "locus", "library", "inFile", "maxIters", "npf", "numThreads")

	if natFile is None or germlineV is None:
		print __doc__
		sys.exit(0)

	#load native sequences
	natives       =  load_fastas(natFile)
	natives_list  =  natives.keys()
	num_nats       =  len(natives_list) + 1


	converged = False #keep track of convergence when running locally


	#load germline sequence
	if not os.path.isfile(library):
		if locus in dict_vgerm_db.keys():
			library = dict_vgerm_db[locus]
		else:
			print "Can't find custom V gene library file!"
			sys.exit(1)

	if inFile != selectedFile:
		correct_V_only = False # Designation of custom starting sequences overrides the -a flag to prevent errors.
		                       # Please do your own filtering if you are using this option.

	germ_dict = load_fastas(library)
	try:
		germ_seq = germ_dict[germlineV]
	except:
		print "Specified germline gene (%s) is not present in the %s library!\n" % (germlineV, library)
		sys.exit(1)
		

	main()

