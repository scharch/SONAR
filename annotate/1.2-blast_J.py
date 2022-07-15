#!/usr/bin/env python3

"""
1.2-blast-J.py

This script parses the BLAST output from 1.1-blast-V_assignment.py. For reads
      for which a V assignment was successfully made, the section of the read
      3' to the V gene is extracted and sent to the cluster for BLAST assignment
      of the J gene. Will also try to assign the D gene if relevant and the
      constant region class.

Usage: 1.2-blast-J.py [ options ]

Options:
    --jlib LIB          Fasta file containing germline J gene sequences. Required only
                            if a custom V gene library was specificied when running
                            to 1.1-blast-V_assignment.py, otherwise the program will use
                            the default libraries matching the V gene locus.
    --dlib LIB          Optional fasta file containing germline D gene sequences, for
                            custom libraries.
    --clib LIB          Optional fasta file containing CH1 gene sequences, for custom
                            libraries.
    --threads 1         Number of threads to use when running locally. [default: 1]
    --cluster           Flag to indicate that blast jobs should be submitted to the
                            SGE cluster. Throws an error if presence of a cluster was
                            not indicated during setup. [default: False]
    --noD               Flag to indicate that blast jobs should be submitted for a
                            D gene library. [default: False]
    --noC               Flag to indicate that no blast jobs should be submitted for a
                            constant region  gene library. [default: False]
    --runFinalize       Flag to call 1.3 script when finished. Additional options to that
                            script (including the ability to automatically call 1.4 and
                            1.5) are listed below. This script will not check the validity
                            of options passed downstream, so user beware. [default: False]

Options for other annotation scripts (see those help messages for details):
    --jmotif TGGGG
    --nterm OPT
    --noclean
    --noFallBack
    --runClustering
    --file FILE
    --min1 1
    --min2 3
    --id .99
    --maxgaps 0
    --runCellStatistics
    --rearrangements rearrangements.tsv
    --save OPT

Created by Zhenhai Zhang on 2011-04-14.
Modified by Chaim A Schramm 2013-07-03 to include j assignment.
Modified by Chaim A Schramm 2014-03-25 to not swamp RAM when processing Illumina
    data.
Edited and commented for publication by Chaim A Schramm on 2015-02-09.
Edited to add queue monitoring by CAS 2015-08-03.
Added local option with threading CAS 2015-11-13.
Edited to use Py3 and DocOpt by CAS 2018-08-22
Updated how Module 1 scripts chain together by CA Schramm 2019-04-01.
Added default constant region DBs for K and L by CAS 2020-01-02.
Updated default database look ups to include `--species` from 1.1 by CAS 2020-02-06.

Copyright (c) 2011-2020 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""

import sys, os, time
from docopt import docopt
from multiprocessing import Pool
from functools import partial

try:
	from SONAR.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/annotate")
	sys.path.append(find_SONAR[0])
	from SONAR.annotate import *


def main():

	if not glob.glob("%s/%s_*.fasta" % (prj_tree.vgene, prj_name)):
		sys.exit("No vBlast output found!\n")

	print( "curating 5'end and strand...." )

	# cut nucleotide sequences from 5'end alignment to germline
	total, good, f_ind = 0, 0, 1
	dict_germ_count	= dict()

	topHandle = open("%s/%s_vgerm_tophit.txt" %(prj_tree.tables, prj_name), "w")
	writer	  = csv.writer(topHandle, delimiter = sep, dialect='unix', quoting=csv.QUOTE_NONE)
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

		print( "%d done, %d good..." %(total, good) )


	f_ind -= 1 #had to go 1 extra to break while loop, now reset to actual number of files
	topHandle.close()


	#print log message
	handle = open("%s/1.2.log" % prj_tree.logs, "w")
	handle.write("total: %d; good: %d\n" %(total, good))
	handle.close()


	#print statistics
	handle = open("%s/%s_vgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
	writer	= csv.writer(handle, delimiter = sep, dialect='unix', quoting=csv.QUOTE_NONE)
	keys	= sorted(dict_germ_count.keys())
	writer.writerow(["gene", "count", "percent"])
	for key in keys:
		aline = [ key, dict_germ_count[key], "%4.2f" % (dict_germ_count[key] / float(total) * 100) ]
		writer.writerow(aline)
	handle.close()


	#run BLAST
	if arguments['--cluster']:

		# write pbs files and auto submit shell script
		if not arguments['--noC']:
			mode = "-db"
			if not os.path.isfile(arguments['--clib'] + ".nhr"):
				mode = "-subject"
			command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s -perc_identity 100" % ( "%03d", CMD_BLAST % (blast_cmd, mode, arguments['--clib'],
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_C_$NUM.txt"   % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
			pbs = open("%s/cblast.sh"%prj_tree.jgene, 'w')
			pbs.write( PBS_STRING%("cBlast-%s"%prj_name, "2G", "1:00:00", "%s 2> %s/%s_C_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
			pbs.close()
			os.system("%s -t 1-%d %s/cblast.sh"%(qsub,f_ind,prj_tree.jgene))

			check = "%s/utilities/checkClusterBlast.py --gene c --big %d --check %s/cmonitor.sh --rehold jMonitor-%s" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene, prj_name)
			monitor = open("%s/cmonitor.sh"%prj_tree.jgene, 'w')
			monitor.write( PBS_STRING%("cMonitor-%s"%prj_name, "2G", "0:30:00", "#$ -hold_jid cBlast-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, check, prj_tree.logs)))
			monitor.close()
			os.system( "%s %s/cmonitor.sh"%(qsub,prj_tree.jgene) )

		if not arguments['--noD']:
			mode = "-db"
			if not os.path.isfile(arguments['--dlib'] + ".nhr"):
				mode = "-subject"
			command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLAST % (blast_cmd, mode, arguments['--dlib'],
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_D_$NUM.txt"   % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
			pbs = open("%s/dblast.sh"%prj_tree.jgene, 'w')
			pbs.write( PBS_STRING%("dBlast-%s"%prj_name, "2G", "1:00:00", "%s 2> %s/%s_D_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
			pbs.close()
			os.system("%s -t 1-%d %s/dblast.sh"%(qsub,f_ind,prj_tree.jgene))

			check = "%s/utilities/checkClusterBlast.py --gene d --big %d --check %s/dmonitor.sh --rehold jMonitor-%s" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene, prj_name)
			monitor = open("%s/dmonitor.sh"%prj_tree.jgene, 'w')
			monitor.write( PBS_STRING%("dMonitor-%s"%prj_name, "2G", "2:00:00", "#$ -hold_jid dBlast-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, check, prj_tree.logs)))
			monitor.close()
			os.system( "%s %s/dmonitor.sh"%(qsub,prj_tree.jgene) )

		#now basic J (do last so the holds work properly -at least for the first round)
		mode = "-db"
		if not os.path.isfile(arguments['--jlib'] + ".nhr"):
			mode = "-subject"
		command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLAST % (blast_cmd, mode, arguments['--jlib'],
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_$NUM.txt"	 % (prj_tree.jgene, prj_name), J_BLAST_WORD_SIZE) )
		pbs = open("%s/jblast.sh"%prj_tree.jgene, 'w')
		pbs.write( PBS_STRING%("jBlast-%s"%prj_name, "2G", "2:00:00", "%s 2> %s/%s_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
		pbs.close()
		os.system("%s -t 1-%d %s/jblast.sh"%(qsub, f_ind,prj_tree.jgene))

		check = "%s/utilities/checkClusterBlast.py --gene j --big %d --check %s/jmonitor.sh" % (SCRIPT_FOLDER, f_ind, prj_tree.jgene)
		if arguments['--runFinalize']:
			check += " --after '%s/annotate/1.3-finalize_assignments.py" % SCRIPT_FOLDER
			for opt in [ '--jmotif', '--nterm', '--file', '--min1', '--min2',
				     '--id', '--maxgaps', '--rearrangements', '--save', '--threads']:
				if arguments[opt] is not None:
					check += " %s %s" % (opt, arguments[opt])
			for flag in ['--cluster', '--noclean', '--noFallBack',
				 '--runClustering', '--runCellStatistics']:
				if arguments[flag]:
					check += " %s" % flag
			check += "'"

		monitor = open("%s/jmonitor.sh"%prj_tree.jgene, 'w')
		monitor.write( PBS_STRING%("jMonitor-%s"%prj_name, "2G", "0:30:00", "#$ -hold_jid jBlast-%s,cMonitor-%s,dMonitor-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, prj_name, prj_name, check, prj_tree.logs))) #wait for C and D to finish before calling 1.3 (if relevant)
		monitor.close()
		os.system( "%s %s/jmonitor.sh"%(qsub,prj_tree.jgene) )

	else:

		#run locally
		partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=arguments['--jlib'], outbase="%s/%s_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE, hits=3)
		blast_pool = Pool(arguments['--threads'])
		blast_pool.map(partial_blast, range(1,f_ind+1))
		blast_pool.close()
		blast_pool.join()

		if not arguments['--noC']:
			partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=arguments['--clib'], outbase="%s/%s_C_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE, hits=3, constant=True)
			blast_pool = Pool(arguments['--threads'])
			blast_pool.map(partial_blast, range(1,f_ind+1))
			blast_pool.close()
			blast_pool.join()

		if not arguments['--noD']:
			partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(prj_tree.jgene, prj_name), db=arguments['--dlib'], outbase="%s/%s_D_%%03d.txt"%(prj_tree.jgene, prj_name), wordSize=J_BLAST_WORD_SIZE)
			blast_pool = Pool(arguments['--threads'])
			blast_pool.map(partial_blast, range(1,f_ind+1))
			blast_pool.close()
			blast_pool.join()

		if arguments['--runFinalize']:
			cmd = "%s/annotate/1.3-finalize_assignments.py" % SCRIPT_FOLDER
			for opt in [ '--jmotif', '--nterm', '--file', '--min1', '--min2',
			             '--id', '--maxgaps', '--rearrangements', '--save', '--threads']:
				if arguments[opt] is not None:
					cmd += " %s '%s'" % (opt, arguments[opt])
			for flag in ['--cluster', '--noclean', '--noFallBack',
		             	 '--runClustering', '--runCellStatistics']:
				if arguments[flag]:
					cmd += " %s" % flag

			print( "Calling 1.3 with command line: %s" % cmd )
			os.system( cmd )




if __name__ == '__main__':

	arguments = docopt(__doc__)
	arguments['--threads'] = int(arguments['--threads'])

	if arguments['--cluster']:
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

	#load saved locus information
	prj_tree	= ProjectFolders(os.getcwd())
	prj_name	= fullpath2last_folder(prj_tree.home)
	with open( "%s/gene_locus.txt" % prj_tree.internal) as handle:
		species = handle.readline().strip()
		locus   = handle.readline().strip()
		vlib    = handle.readline().strip()

	#check germline libraries
	if arguments['--jlib'] is not None:
		if not os.path.isfile(arguments['--jlib']):
			print( "Can't find custom J gene library file!" )
			sys.exit(1)
	elif locus == "C":
		print( "Error: Custom V gene library was used; please specify matching J gene library" )
		sys.exit(1)
	else:
		arguments['--jlib'] = eval( SUPPORTED_SPECIES[species] + "_J" + locus + "_DB" )

	if not arguments['--noD']:
		if arguments['--dlib'] is not None:
			if not os.path.isfile(arguments['--dlib']):
				print( "Can't find custom D gene library file!" )
				sys.exit(1)
		elif locus == "C":
			print( "Error: Custom V gene library was used; please specify matching D gene library or use the --noD flag" )
			sys.exit(1)
		elif "H" in locus: #include HKL
			arguments['--dlib'] = eval( SUPPORTED_SPECIES[species] + "_DH_DB" )
		else:
			arguments['--noD'] = True #light chains only

	if not arguments['--noC']:
		if arguments['--clib'] is not None:
			if not os.path.isfile(arguments['--clib']):
				print( "Can't find custom constant region gene library file!" )
				sys.exit(1)
		elif locus == "C":
			print( "Error: Custom V gene library was used; please specify matching constant region gene library or use the --noC flag" )
			sys.exit(1)
		else:
			arguments['--clib'] = eval( SUPPORTED_SPECIES[species] + "_C" + locus + "_DB" )


	#log command line
	logCmdLine(sys.argv)

	# save J/D/C library locations for next step
	# overwrite the file to avoid errors in 1.3 in case we are experimenting with multiple germline databases
	with open( "%s/gene_locus.txt" % prj_tree.internal, 'w') as handle:
		handle.write("%s\n" % species)
		handle.write("%s\n" % locus)
		handle.write("%s\n" % vlib)
		handle.write("%s\n" % arguments['--jlib'])
		handle.write("%s\n" % arguments['--dlib'])
		handle.write("%s\n" % arguments['--clib'])
		handle.close()


	main()
