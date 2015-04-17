#!/usr/bin/python

"""
2.4-lineage-intradonor.py

This script uses an iterative phylogenetic analysis to find sequences related
      to a set of known antibodies. Preprocessed sequences (including optional
      filtering by assigned germline) are split into groups of 250 and used to
      together with known sequences to build neighbor-joining trees rooted on
      the germline V gene of the known antibodies. Sequences in the minimum
      sub-tree spanning all of the known sequences are passed forward into the
      next iteration. The algorithm is considered to have converged when 95%
      of the input sequences in a round are in the minimum sub-tree.

This script has an option to run the analysis as a cluster job, in which case
      2,500 sequences will be processed in each group.

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

Usage: 2.4-lineage-intradonor.py -n native.fa -v germline_V
                                 [-locus <H|K|L|C> -lib path/to/library.fa
				  -i custom/input.fa -maxIters 10
				  -nofilter -a -cluster -h -f]

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
    maxIters	Optional maximum number of rounds to conduct before giving up.
                   Default = 10.

    Optional flags:
    nofilter	Do NOT filter NGS sequences for correct germline V gene
                   assignment. Default = OFF (DO filter).
    a		Use all NGS sequences with an assigned V, even those with
                   out-of frame junctions and/or stop codons or without a 
		   successfuly assigned J gene. Default = OFF (use in-frame
		   ORF sequences only).
    cluster	Submit tree-building jobs to the cluster.
    f		Force a restart of the analysis, even if there are files from
                   a previous run in the working directory.

Created by Zhenhai Zhang on 2011-07-12.
Modified to compress into a single script and many updates by 
        Chaim A Schramm 2014-01-09.
Edited and commented for publication by Chaim A Schramm on 2015-04-14.

Copyright (c) 2011-2015 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import time, sys

#assume that we haven't necessarily set the path variables on the cluster
find_ZAP_on_cluster = sys.argv[0].split("zap")
sys.path.append("%szap" % find_ZAP_on_cluster[0])

from cStringIO import StringIO
from zap.lineage import *
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline



def main():

	global npf, converged, maxIters, cluster, force, num_nats, correct_V_only, germlineV, germ_seq, inFile, natFile, locus, library
	currentIter = 0

	# master loop
	while not converged:
		
		# parse tree files and get all reads clustered with native antibodies
		# start total at 1 so good/total < .95 when we start a new analysis
		tree_files, good, total, retained_reads = sorted(glob.glob("%s/NJ*.tree" %prj_tree.clustal)), 0.0, 1.0, []
								 
		if force:
			tree_files = [] #no need to process the files since we will be starting over
			force = False #turn it off so we don't get stuck in an infinite loop of restarts

	        #this gets skipped in the first round
		for tf in tree_files:
			
		        #read in the tree as a string
			tree_string = open(tf, "rU").read().strip()
			
			#get rid of possible negative branch lengths and then parse
			tree_string, num = re.subn(":-\d+\.\d+",":0", tree_string)
			tree_string, num = re.subn("\n","", tree_string)
			tree = Phylo.read(StringIO(tree_string), "newick")

			
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
			#print "Found %d reads in subtree. Total saved: %d / %d" % ( len(all_leaves)-num_nats, good, total-1)


		# processed all trees from last round, now do a sanity check
		if good == 0 and len(tree_files) > 0:
			sys.exit( "NO positive sequences found --stopped!" )

		
		# Are we starting a new run? If not, limit memory usage by only reading in the retained reads
		read_dict = dict()
		if len(tree_files) == 0:
			print "%s - Starting a new analysis from scratch..." % time.strftime("%H:%M:%S")
			if correct_V_only:
				read_dict = load_fastas_with_Vgene( inFile, germlineV )
			else:
				read_dict = load_fastas( inFile )
		else:
			read_dict = load_seqs_in_dict( inFile, retained_reads )
		#randomize the order
		shuffled_reads = read_dict.values()
		random.shuffle(shuffled_reads)
			
		# Check for convergence before starting a new round
		if good/total < 0.95:
			
			if currentIter >= maxIters:
				sys.exit( "Maximum number of iterations reached without convergence. Current round: %d reads, %5.2f%% of input" % (good, good/total) )
			else:
				if len(tree_files) > 0: 
					#it's a silly message to print the first time through
					print "%s - Finished processing round %d: %d reads, %5.2f%% of input" % (time.strftime("%H:%M:%S"), currentIter, good, 100*good/total)
				currentIter +=1


			#do some cleanup
			oldFiles = glob.glob("%s/*" % prj_tree.last)
			for old in oldFiles:
				os.remove(old)
			lastRound = glob.glob("%s/NJ*" % prj_tree.clustal)
			for infile in lastRound:
				os.rename( infile, "%s/%s" % (prj_tree.last, os.path.basename(infile)) )
					
			#open initial output
			f_ind = 1
			outFasta = open("%s/NJ%05d.fa" % (prj_tree.clustal, f_ind), "w")
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
					if not cluster:
						run_muscle            = MuscleCommandline( input="%s/NJ%05d.fa" % (prj_tree.clustal, f_ind), out="%s/NJ%05d.aln" % (prj_tree.clustal, f_ind) )
						run_muscle.tree1      = "%s/NJ%05d.tree" % (prj_tree.clustal, f_ind)
						run_muscle.cluster1   = "neighborjoining"
						run_muscle.maxiters   = 1
						thisVarHidesTheOutput = run_muscle()
						if f_ind % 25 == 0: print "%s - Finished %dth file of %d in this round..." % (time.strftime("%H:%M:%S"), f_ind, len(shuffled_reads)/npf + 1)
					f_ind += 1
					outFasta = open("%s/NJ%05d.fa" % (prj_tree.clustal, f_ind), "w")
					currentSize = 0

			#Any leftovers?
			if currentSize > 0:
				outFasta.write( ">%s\n%s\n" % (germ_seq.id, germ_seq.seq) )
				for n in natives.values():
					outFasta.write( ">%s\n%s\n" % (n.id, n.seq) )
				outFasta.close()
			
				if not cluster:
					run_muscle = MuscleCommandline( input="%s/NJ%05d.fa" % (prj_tree.clustal, f_ind), out="%s/NJ%05d.aln" % (prj_tree.clustal, f_ind) )
					run_muscle.tree1      = "%s/NJ%05d.tree" % (prj_tree.clustal, f_ind)
					run_muscle.cluster1   = "neighborjoining"
					run_muscle.maxiters   = 1
					thisVarHidesTheOutput = run_muscle()
			else:
				outFasta.close()
				f_ind -= 1 #don't submit an empty file to the cluster

				
			# At this point, if we are running on the cluster, submit this round and quit loop
			if cluster:
				# write pbs files and auto submit shell script
				command = "NUM=`printf \"%%05d\" $SGE_TASK_ID`\n%s -in %s/NJ$NUM.fa -out %s/NJ$NUM.aln -cluster1 neighborjoining -maxiters 1 -tree1 %s/NJ$NUM.tree" % \
				    (cluster_muscle, prj_tree.clustal, prj_tree.clustal, prj_tree.clustal)
				pbs = open("%s/intradonor.sh" % prj_tree.clustal, 'w')
				pbs.write( PBS_STRING%(f_ind, "%s-intradonor"%prj_name, "500M", "10:00:00", "%s > %s/NJ$NUM.out 2> %s/NJ$NUM.err"%(command, prj_tree.clustal, prj_tree.clustal)) )
				pbs.close()

				# write a second PBS job to restart this script once muscle has finished
				next_round = open("%s/nextround.sh" % prj_tree.clustal, 'w')
				next_round.write( "\
#!/bin/bash\n\
#$ -hold_jid %s-intradonor\t\t# run once previous round has completed\n\
#$ -N restart-intradonor\t\t# job name\n\
#$ -l mem=1G,time=1:00:00\t\t# resource requests\n\
#$ -cwd\t\t\t\t\t# use current directory\n\
#$ -o %s/restart-intradonor.out\t#output\n\
#$ -e %s/restart-intradonor.err\t#error\n\
%s/lineage/2.4-intradonor.py -n %s -v %s -locus %s -lib %s -i %s -maxIters %d -cluster\n" % 
						  (prj_name, prj_tree.clustal, prj_tree.clustal, ZAP_FOLDER, natFile, germlineV, locus, library, inFile, maxIters-1) )
				next_round.close()

				os.system("qsub %s/intradonor.sh" % prj_tree.clustal)
				os.system("qsub %s/nextround.sh"  % prj_tree.clustal)
				break #exits "while not converged" loop

		
		# If we have converged, write final output
		else:
			out_list = open("%s/%s_intradonor_positives.txt" % (prj_tree.tables, prj_name), "w")
			out_nt   = open("%s/%s_intradonor_positives.fa"  % (prj_tree.nt, prj_name),     "w")
			out_aa   = open("%s/%s_intradonor_positives.fa"  % (prj_tree.aa, prj_name),     "w")
			
			for r in read_dict.values():
				out_list.write( ">%s\n%s\n" % (r.id, r.seq) )
				out_nt.write( ">%s\n%s\n" % (r.id, r.seq) )
				out_aa.write( ">%s\n%s\n" % (r.id, r.seq.translate()) )
				
			out_list.close()
			out_nt.close()
			out_aa.close()
			
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
	selectedFile = "%s/%s_goodVJ.fa" % (prj_tree.nt, prj_name)
	flag = [x for x in ["a", "-a", "--a", "all", "-all", "--all"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		selectedFile = "%s/%s_allV.fa" % (prj_tree.nt, prj_name)

	#check cluster parameter
	cluster = False
	npf = 250
	flag = [x for x in ["c", "-c", "--c", "cluster", "-cluster", "--cluster"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		cluster = True
		npf = 2500

	#get the parameters from the command line
	dict_args = processParas(sys.argv, n="natFile", v="germlineV", locus="locus", lib="library", i="inFile", maxIters="maxIters")
	defaults = dict(locus="H", library="", inFile=selectedFile, maxIters=10)
	natFile, germlineV, locus, library, inFile, maxIters = getParasWithDefaults(dict_args, defaults, "natFile", "germlineV", "locus", "library", "inFile", "maxIters")

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

