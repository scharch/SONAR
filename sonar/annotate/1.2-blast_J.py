#!/usr/bin/env python

"""
1.2-blast-J.py

This script parses the BLAST output from 1.1-blast-V_assignment.py. For reads
      for which a V assignment was successfully made, the section of the read
      3' to the V gene is extracted and sent to the cluster for BLAST assignment
      of the J gene. Will also try to assign the D gene if relevant and the 
      constant region class. 

Usage: 1.2-blast-J.py -lib  path/to/j-library.fa
                      -dlib path/to/d-library.fa
		      -clib path/to/c-library.fa
		      -threads 1 -cluster
		      -noD -noC
		      -callFinal -h

    All parameters are optional. Invoke with -h or --help to print this
        documentation.

    lib 	fasta file containing germline J gene sequences. Required only
                     if "-locus C" was specificied to 1.1-blast-V_assignment.py;
		     otherwise the program will use the default libraries.
    dlib 	Optional fasta file containing germline D gene sequences, for 
                     custom libraries.
    clib 	Optional fasta file containing CH1 gene sequences, for custom
                     libraries.
    threads     Number of threads to use when running locally. Ignored if 
                   -cluster is specified. Default = 1.
    cluster     Flag to indicate that blast jobs should be submitted to the
                   SGE cluster. Throws an error if presence of a cluster was
		   not indicated during setup. Default = run locally.
    noD         Flag to indicate that no blast jobs should be submitted for a
                   D gene library. Default = False (do D gene blast) unless a 
		   light chain library is specified.
    noC         Flag to indicate that no blast jobs should be submitted for a
                   constant region  gene library. Default = False (do constant
		   region blast) unless a light chain library is specified.
    callFinal   Optional flag to call 1.3-finalize_assignments.py when done.
                     Default = False.

Created by Zhenhai Zhang on 2011-04-14.
Modified by Chaim A Schramm 2013-07-03 to include j assignment.
Modified by Chaim A Schramm 2014-03-25 to not swamp RAM when processing Illumina
    data.
Edited and commented for publication by Chaim A Schramm on 2015-02-09.
Edited to add queue monitoring by CAS 2015-08-03.
Added local option with threading CAS 2015-11-13.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""

import sys, os, time
from multiprocessing import Pool
from functools import partial

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *


def main():
	
        if not glob.glob("%s/%s_*.fasta" % (prj_tree.vgene, prj_name)):
                sys.exit("No vBlast output found!\n")
        
	print "curating 5'end and strand...."

	# cut nucleotide sequences from 5'end alignment to germline
	total, good, f_ind = 0, 0, 1
	dict_germ_count	= dict()
	
        topHandle = open("%s/%s_vgerm_tophit.txt" %(prj_tree.tables, prj_name), "w")
	writer    = csv.writer(topHandle, delimiter = sep)
	writer.writerow(PARSED_BLAST_HEADER)
	
	while os.path.isfile("%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind)):

		split_fasta = "%s/%s_%03d.fasta" %(prj_tree.jgene, prj_name, f_ind)
		fasta_handle = open(split_fasta, "w")

		# parse blast output
		dict_germ_aln, dict_other_germs, dict_germ_count = get_top_hits( "%s/%s_%03d.txt" % (prj_tree.vgene, prj_name, f_ind), topHitWriter=writer, dict_germ_count=dict_germ_count )
	
		# process each sequence
		for entry in SeqIO.parse("%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind), "fasta"):

			total += 1
			if entry.id in dict_germ_aln:
				if dict_germ_aln[entry.id].strand == "plus":
					entry.seq = entry.seq[ dict_germ_aln[entry.id].qend : ]
				else:
					entry.seq = entry.seq[ : dict_germ_aln[entry.id].qstart -1 ]
					entry.seq = entry.reverse_complement().seq

				if len(entry.seq) > 30: #can probably be 50...
					fasta_handle.write(">%s\n%s\n" % (entry.id,entry.seq))
					good += 1

		fasta_handle.close()
		f_ind += 1
		
		print "%d done, %d good..." %(total, good)


	f_ind -= 1 #had to go 1 extra to break while loop, now reset to actual number of files
        topHandle.close()
        

	#print log message
	handle = open("%s/1.2.log" % prj_tree.logs, "w")
	handle.write("total: %d; good: %d\n" %(total, good))
	handle.close()
	

        #print statistics
	handle = open("%s/%s_vgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
	writer 	= csv.writer(handle, delimiter = sep)
	keys 	= sorted(dict_germ_count.keys())
	writer.writerow(["gene", "count", "percent"])
	for key in keys:
		aline = [ key, dict_germ_count[key], "%4.2f" % (dict_germ_count[key] / float(total) * 100) ]
		writer.writerow(aline)
	handle.close()


	#run BLAST
	if useCluster:

		# write pbs files and auto submit shell script
		if os.path.isfile(const_lib):
			command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s -perc_identity 100" % ( "%03d", CMD_BLAST % (cluster_blast, const_lib, 
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_C_$NUM.txt"   % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
			pbs = open("%s/cblast.sh"%prj_tree.jgene, 'w')
			pbs.write( PBS_STRING%("cBlast-%s"%prj_name, "2G", "1:00:00", "%s 2> %s/%s_C_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
			pbs.close()
			os.system("%s -t 1-%d %s/cblast.sh"%(qsub,f_ind,prj_tree.jgene))

			check = "%s/utilities/checkClusterBlast.py -gene c -big %d -check %s/cmonitor.sh -rehold jMonitor%s" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene, prj_name)
			monitor = open("%s/cmonitor.sh"%prj_tree.jgene, 'w')
			monitor.write( PBS_STRING%("cMonitor-%s"%prj_name, "2G", "0:30:00", "#$ -hold_jid cBlast-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, check, prj_tree.logs)))
			monitor.close()
			os.system( "%s %s/cmonitor.sh"%(qsub,prj_tree.jgene) )

		if os.path.isfile(dlib):
			command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLAST % (cluster_blast, dlib, 
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_D_$NUM.txt"   % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
			pbs = open("%s/dblast.sh"%prj_tree.jgene, 'w')
			pbs.write( PBS_STRING%("dBlast-%s"%prj_name, "2G", "1:00:00", "%s 2> %s/%s_D_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
			pbs.close()
			os.system("%s -t 1-%d %s/dblast.sh"%(qsub,f_ind,prj_tree.jgene))

			check = "%s/utilities/checkClusterBlast.py -gene d -big %d -check %s/dmonitor.sh -rehold jMonitor%s" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene, prj_name)
			monitor = open("%s/dmonitor.sh"%prj_tree.jgene, 'w')
			monitor.write( PBS_STRING%("dMonitor-%s"%prj_name, "2G", "2:00:00", "#$ -hold_jid dBlast-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, check, prj_tree.logs)))
			monitor.close()
			os.system( "%s %s/dmonitor.sh"%(qsub,prj_tree.jgene) )

		#now basic J (do last so the holds work properly -at least for the first round)
		command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLAST % (cluster_blast, library, 
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_$NUM.txt"   % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
		pbs = open("%s/jblast.sh"%prj_tree.jgene, 'w')
		pbs.write( PBS_STRING%("jBlast-%s"%prj_name, "2G", "2:00:00", "%s 2> %s/%s_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
		pbs.close()
		os.system("%s -t 1-%d %s/jblast.sh"%(qsub, f_ind,prj_tree.jgene))

		check = "%s/utilities/checkClusterBlast.py -gene j -big %d -check %s/jmonitor.sh" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene)
		if callF:
			check += " -after %s/annotate/1.3-finalize_assignments.py" % SCRIPT_FOLDER
		monitor = open("%s/jmonitor.sh"%prj_tree.jgene, 'w')
		monitor.write( PBS_STRING%("jMonitor-%s"%prj_name, "2G", "0:30:00", "#$ -hold_jid jBlast-%s,cMonitor-%s,dMonitor-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, prj_name, prj_name, check, prj_tree.logs))) #wait for C and D to finish before calling 1.3 (if relevant)
		monitor.close()
		os.system( "%s %s/jmonitor.sh"%(qsub,prj_tree.jgene) )

	else:

		#run locally
                partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=library, outbase="%s/%s_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE, hits=3)
                blast_pool = Pool(numThreads)
                blast_pool.map(partial_blast, range(1,f_ind+1))
                blast_pool.close()
                blast_pool.join()

		if os.path.isfile(const_lib):
			partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=const_lib, outbase="%s/%s_C_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE, hits=3, constant=True)
			blast_pool = Pool(numThreads)
			blast_pool.map(partial_blast, range(1,f_ind+1))
			blast_pool.close()
			blast_pool.join()

		if os.path.isfile(dlib):
			partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=dlib, outbase="%s/%s_D_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE)
			blast_pool = Pool(numThreads)
			blast_pool.map(partial_blast, range(1,f_ind+1))
			blast_pool.close()
			blast_pool.join()

		if callF:
			os.system( "%s/annotate/1.3-finalize_assignments.py" % SCRIPT_FOLDER )




if __name__ == '__main__':
	
	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	#log command line
	logCmdLine(sys.argv)

	#check cluster usage
	useCluster = False
	if q("-cluster"):
		sys.argv.remove("-cluster")
		useCluster = True
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

        #check if call 1.3
        callF = False
        if q("-callFinal"):
                sys.argv.remove("-callFinal")
                callF = True

        #check if blast D
        blastD = True
        if q("-noD"):
                sys.argv.remove("-noD")
                blastD = False

        #check if blast C
        blastC = True
        if q("-noC"):
                sys.argv.remove("-noC")
                blastC = False

	#get parameters from input
	dict_args = processParas(sys.argv, lib="library", dlib="dlib", clib="const_lib", threads = "numThreads")
	library, dlib, const_lib, numThreads = getParasWithDefaults(dict_args, dict(library="", dlib="", const_lib="", numThreads=1), "library", "dlib", "const_lib", "numThreads")

	
	prj_tree        = ProjectFolders(os.getcwd())
	prj_name        = fullpath2last_folder(prj_tree.home)


	#load saved locus information
	handle = open( "%s/gene_locus.txt" % prj_tree.internal)
	locus = handle.readline().strip()
	handle = open( "%s/gene_locus.txt" % prj_tree.internal, 'a+')
	# we'll keep custom libraries even for a default locus (maybe someone wants to use an updated set of D alleles?)
	if not os.path.isfile(library):
		if locus in dict_jgerm_db.keys():
			library = dict_jgerm_db[locus]
		else:
			print "Can't find custom J gene library file!"
			sys.exit(1)
	if locus == "H":
		if blastD and not os.path.isfile(dlib)     : dlib      = DH_DB
		if blastC and not os.path.isfile(const_lib): const_lib = CH_DB

	# save J/D/C library locations for next step
	handle.write("%s\n" % library)
	handle.write("%s\n" % dlib)
	handle.write("%s\n" % const_lib)
	handle.close()

	main()

