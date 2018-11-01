#!/usr/bin/env python3

"""
1.1-blast-V.py

This script looks for raw NGS data files in the current folder and parses them
      into manageable chunks for searching with BLAST on the cluster. Each raw
      read that meets the length restrictions is assigned a unique 8 digit
      serial number and submitted to the cluster for BLASTing against the 
      desired V gene germline library.
      Distinctions between raw data files will not be maintained, although they
      can be reconstructed using the id table. The name of the current folder
      is assigned as the project name, which will be used to identify all
      output files created.

Usage: 1.1-blast-V.py [ options ] [ --locus H | --lib path/to/library.fa ] [ --fasta file1.fa ]... [ --derep ] [ --cluster | [--npf 10000 --threads 1] ] 

Options:
    --locus H		 H: heavy chain / K: kappa chain / L: lambda chain / 
			    KL: kappa OR lambda / HKL: any. [default: H]
    --lib LIB            Location of file containing custom library (e.g. for use with
			    non-human genes). Mutually exclusive with  --locus.
    --fasta <input.fa>	 File(s) containing the input reads to process. May be specified
			    multiple times. By default, uses all FASTA/FASTQ files (those
			    with extensions of .fa, .fas, .fst, .fasta, .fna, .fq, or
			    .fastq) in the root project directory.
    --derep              Flag to dereplicate input sequences using vsearch prior to 
                            processing for blast. The number of reads supporting each
                            sequence will be saved and output as "duplicate_count."
                            Results will also be available in 
                            work/internal/derepAllRawSeqs.uc [default: False]
    --cluster            Flag to indicate that blast jobs should be submitted to the
                            SGE cluster. Throws an error if presence of a cluster was
                            not indicated during setup. [default: False]
    --npf <10000>	 Number of sequences in each split file when running locally. 
			    (Resource requests for a cluster are calibrated to groups
			    of 50K sequences, and cannot be changed.) [default: 10000]
    --threads <1>        Number of threads to use when running locally. [default: 1]
    --callJ              Flag to call 1.2-blast_J.py when done. [default: False]
    --jArgs <options>	 Optional arguments to be provided to 1.2-blast_j.py. If provided,
			    forces callJ flag to True.
    --minl <300>         Minimum length for read filtering (inclusive). [default: 300]
    --maxl <600>         Maximum length for read filtering (inclusive). [default: 600]
    -f                   Force new analysis and overwrite existing working directories.
                            [default: False]
    --qual		 CURRENTLY DEPRECATED! Use PHRED scores for QC. [default: False]

Created by Zhenhai Zhang on 2011-04-12.
Edited and commented for publication by Chaim A Schramm on 2014-12-22.
Edited to add queue monitoring by CAS 2015-07-30.
Added local option with threading CAS 2015-11-13.
Edited to use Py3 and DocOpt by CAS 2018-08-22.
Added derep option by CAS 2018-10-19.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
			 Institutes of Health, USA. All rights reserved.

"""


import sys
import os
import time
from docopt import docopt
from multiprocessing import Pool
from functools import partial

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *



def getSeqsWithFileName(inFile):
	#helper function to add source file info to a seqID
	#  to preserve info through derep
	filetype = "fasta"
	if re.search("\.(fq|fastq)$", inFile) is not None:
		filetype = "fastq"
		with open(inFile, "rU") as handle:
			for entry in SeqIO.parse(handle, filetype):
				entry.id += ";file=%s"%inFile
				yield entry



def main():

	#if no input files were specified, glob up everything
	if len(arguments['--fasta'])==0: arguments['--fasta'] = glob.glob("*.fa") + glob.glob("*.fas") + glob.glob("*.fst") + glob.glob("*.fasta") + glob.glob("*.fna") + glob.glob("*.fq") + glob.glob("*.fastq")

	#dereplicate?
	if arguments['--derep']:
		with open( "%s/tempForDerep.fa"%folder_tree.internal, "w" ) as tempFile:
			for myFile in arguments['--fasta']:
				SeqIO.write( getSeqsWithFileName(myFile), tempFile, "fasta" )

		subprocess.call( [ vsearch, "-derep_fulllength", "%s/tempForDerep.fa"%folder_tree.internal,
				   "-output", "%s/temp_derep.fa"%folder_tree.internal,
				   "-uc", "derepAllRawSeqs.uc",
				   "-sizein", "-sizeout" ] )

		arguments['--fasta'] = [ "%s/temp_derep.fa"%folder_tree.internal ]
		os.remove( "%s/tempForDerep.fa"%folder_tree.internal )
				

	#initiate counters
	total, total_good, f_ind = 0, 0, 1

	
	# open output files
	fasta	  = open("%s/%s_%03d.fasta"  % (folder_tree.vgene,  prj_name, f_ind), 'w')
	id_handle = open("%s/id_lookup.txt"  %	folder_tree.internal,		      'w')
	id_map	  = csv.writer(id_handle, delimiter=sep)


	#if we decide to use quals for something, can add this block back in
	'''
	if arguments['--qual']:
		qual =	     open("%s/%s_%03d.qual"   % (folder_tree.vgene,    prj_name, f_ind), 'w')
		qual_generator	= generate_quals_folder(folder_tree.home)
	'''


	#iterate over input and split for blast
	for eachFile in arguments['--fasta']:
		for sequence in generate_read_fasta( eachFile ):

			sourceFile = eachFile
			fromDerep  = re.search(";file=([^;\s]+)", sequence.id)
			if fromDerep:
				sourceFile = fromDerep.group(1)

			dup_count = 1
			checkSize = re.search(";size=(\d+)", sequence.id)
			if checkSize:
				dup_count = checkSize.group(1)
				
			total += 1
			id_map.writerow([ "%08d"%total, sourceFile, sequence.id, dup_count, len(sequence.seq)])

			if arguments['--minl'] <= len(sequence.seq) <= arguments['--maxl']:
				total_good += 1
				fasta.write(">%08d\n%s\n" % (total, sequence.seq))

				#uncomment to re-implement quals
				'''
				if arguments['--qual']:
				if re.search("\.(fq|fastq)$", file_name) is None:
					myqual = qual_generator.next()
				qual.write(">%08d\n%s\n" % (total, " ".join(map(str, myqual.qual_list))))
				'''

				if total_good % arguments['--npf'] == 0: 
					#close old output files, open new ones, and print progress message
					fasta.close()
					f_ind += 1
					fasta = open("%s/%s_%03d.fasta" % (folder_tree.vgene, prj_name, f_ind), 'w')
					print( "%d processed, %d good; starting file %s_%03d" %(total, total_good, prj_name, f_ind) )

					'''
					if arguments['--qual']:
					qual.close()
					qual = open("%s/%s_%03d.qual"%(folder_tree.vgene, prj_name, f_ind), 'w')
					'''
				
	print( "TOTAL: %d processed, %d good" %(total, total_good) )
	
	fasta.close()
	id_handle.close()
	'''
	if arguments['--qual']:
		qual.close()
	'''


	#clean up - but leave the uc file, in case we want to track replicates later
	if arguments['--derep']:
		os.remove( "%s/temp_derep.fa"%folder_tree.internal )

		
	#print log message
	handle = open("%s/1-split.log" % folder_tree.logs, "w")
	handle.write("total: %d; good: %d; percentile: %f\n" %(total, total_good, float(total_good)/total * 100))
	handle.close()
	

	#run BLAST
	if arguments['--cluster']:

		# write pbs files and auto submit shell script
		command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLAST % (blast_cmd, arguments['--lib'], 
									      "%s/%s_$NUM.fasta" % (folder_tree.vgene, prj_name),
									      "%s/%s_$NUM.txt"	 % (folder_tree.vgene, prj_name), V_BLAST_WORD_SIZE) )
		pbs = open("%s/vblast.sh"%folder_tree.vgene, 'w')
		pbs.write( PBS_STRING%("vBlast-%s"%prj_name, "2G", "2:00:00", "%s 2> %s/%s_$NUM.err"%(command, folder_tree.vgene, prj_name)) )
		pbs.close()
		os.system( "%s -t 1-%d %s/vblast.sh"%(qsub,f_ind,folder_tree.vgene) )
		
		check = "%s/utilities/checkClusterBlast.py --gene v --big %d --check %s/vmonitor.sh" % (SCRIPT_FOLDER, f_ind, folder_tree.vgene)
		if arguments['--callJ']:
			check += " --after '%s/annotate/1.2-blast_J.py %s'" % (SCRIPT_FOLDER, arguments['--jArgs'])
			monitor = open("%s/vmonitor.sh"%folder_tree.vgene, 'w')
			monitor.write( PBS_STRING%("vMonitor-%s"%prj_name, "2G", "1:00:00", "#$ -hold_jid vBlast-%s\n%s >> %s/qmonitor.log 2>&1"%(prj_name, check, folder_tree.logs)) )
			monitor.close()
			os.system( "%s %s/vmonitor.sh"%(qsub,folder_tree.vgene) )

	else:

		#run locally
		partial_blast = partial( blastProcess, filebase="%s/%s_%%03d.fasta"%(folder_tree.vgene, prj_name), db=arguments['--lib'], outbase="%s/%s_%%03d.txt"%(folder_tree.vgene, prj_name), wordSize= V_BLAST_WORD_SIZE)
		blast_pool = Pool(arguments['--threads'])
		blast_pool.map(partial_blast, range(1,f_ind+1))
		blast_pool.close()
		blast_pool.join()

		if arguments['--callJ']:
			os.system( "%s/annotate/1.2-blast_J.py %s" % (SCRIPT_FOLDER, arguments['--jArgs']) )


	

if __name__ == '__main__':

	#save command line
	cmdLine = sys.argv

	arguments = docopt(__doc__)

	arguments['--threads'] = int(arguments['--threads'])
	arguments['--npf']     = int(arguments['--npf'])
	arguments['--minl']    = int(arguments['--minl'])
	arguments['--maxl']    = int(arguments['--maxl'])
	
	
	if arguments['--cluster']:
		arguments['--npf'] = 50000
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

	#check if call J
	if arguments['--jArgs'] is not None:
		arguments['--callJ'] = True
	elif arguments['--callJ']:
		arguments['--jArgs'] = "" #kludge to prevent errors

	if arguments['--lib'] is not None:
		arguments['--locus'] = 'C'
		if not os.path.isfile(arguments['--lib']):
			print( "Can't find custom V gene library file!" )
			sys.exit(1)                
	else:
		if arguments['--locus'] in dict_vgerm_db.keys():
			arguments['--lib'] = dict_vgerm_db[arguments['--locus']]
		else:
			print("Error: valid options for --locus are H, K, L, KL, and HKL only")
			sys.exit(1)

	# create 1st and 2nd subfolders
	prj_folder  = os.getcwd()
	folder_tree = create_folders( prj_folder, force=arguments['-f'] )
	prj_name    = prj_folder[prj_folder.rindex("/") + 1 :]

	#log command line
	logCmdLine(sys.argv)	

	handle = open( "%s/gene_locus.txt" % folder_tree.internal, "w")
	handle.write( "%s\n%s\n" % (arguments['--locus'], arguments['--lib']) )
	handle.close()

	main()

