#!/usr/bin/env python3

"""
1.3-finalize_assignments.py

This script parses the BLAST output from 1.1-blast-V_assignment.py and 
      1.2-blast_J.py. Sequences with successful assignments are
      output into fasta files and a master table is created summarizing the
      properties of all input sequences.

Usage:  1.3-finalize_assignments.py [ --jmotif "TT[C|T][G|A]G" --nterm truncate --noclean ] [ options ]

Options:
    --jmotif TGGGG    Conserved nucleotide sequence indicating the start of FWR4 on
                          the J gene. Defaults to either TGGGG for heavy chains or
                          TT[C|T][G|A]G for light chains; set manually for custom 
                          light chain libraries or for species with a different motif.
    --nterm OPT       What to do if blast hit does not extend to the N terminus of the
                          germline V gene. Options are:
                              truncate - trim to blast hit;
                              extend   - change trim boundary to correspond to expected N
                                             N terminus of germline or beginning of read if
                                             shorter (NOT recommended, typically results in
                                             bad sequences)
                              germline - replace missing region with the germline V
                                             sequence (useful for FWR1 primers)
                              discard  - only mark as 'good' sequences with full germline
                                             region
                          [default: truncate]
    --noclean         Disable automatic deletion of working files from 1.1 and 1.2 [default: False]
    --noFallBack      Flag to disable fall-back detection of heavy chain isotype based on the
                          first 3 bases of CH1. Useful for species where those may be different
                          than in humans. [default: False]
    --cluster         Flag to indicate that blast jobs should be submitted to the
                          SGE cluster. Throws an error if presence of a cluster was
                          not indicated during setup. [default: False]
    --threads <1>     Number of threads to use when running locally. [default: 1]
    --reenter         Flag used internally to resume after submitting individual jobs
                          to the cluster. [default: False]
    --runClustering   Flag to call 1.4 script when finished. Additional options to that
                          script (including the ability to automatically call 1.5) are listed
                          below. This script will not check the validity of options passed
                          downstream, so user beware. [default: False]

Options for other annotation scripts (see those help messages for details):
    --file FILE 
    --min1 1
    --min2 3
    --id .99
    --maxgaps 0
    --runCellStatistics
    --rearrangements rearrangements.tsv
    --save OPT

Created by Chaim A Schramm on 2013-07-05
Edited and commented for publication by Chaim A Schramm on 2015-02-25.
Edited to add custom J motif option for other species by CAS 2016-05-16.
Edited to add options to handle missing N terminal by CAS 2017-05-08.
Edited to make motif matching case-insensitive and to disable auto-clean-up
        by CAS 2017-07-31.
Edited to distinguish between nonproductive junctions and apparent frameshifts
        within V by CAS 2018-07-23.
Edited to use Py3 and DocOpt by CAS 2018-08-22.
Edited to use AIRR datarep convention by CAS 2018-10-12.
Added pulling through of cell/umi stats from 1.0 by CAS 2019-03-01.
Added fall-back checking of isotype (IgG vs IgM) based on first 3 bases of CH1
        since this is available with current Douek lab primers by CAS 2019-03-15. 
Updated how Module 1 scripts chain together by CA Schramm 2019-04-01.
Moved main processing code to parse_blast.py to allow for parallelization on a
        cluster by CA Schramm 2019-04-01.

Copyright (c) 2011-2019 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""


import sys, os
from docopt import docopt
import airr
from Bio import Seq
from multiprocessing import Pool
from collections import Counter

try:
	from SONAR.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/annotate")
	sys.path.append(find_SONAR[0])
	from SONAR.annotate import *
	
def callParser(chunk):
	cmd = "%s/annotate/parse_blast.py --jmotif '%s' --nterm %s --chunk %03d" % \
						( SCRIPT_FOLDER, arguments['--jmotif'], arguments['--nterm'], chunk )
	if arguments['--noFallBack']: cmd += " --noFallBack"
	os.system( cmd )


def main():

	if not glob.glob("%s/%s_*.fasta" % (prj_tree.jgene, prj_name)):
		sys.exit("No jBlast output found!\n")
		
	maxFiles = len( glob.glob("%s/%s_*.fasta" % (prj_tree.vgene, prj_name)) )
	
	if not arguments['--reenter']:
		print( "curating junction and 3' end..." )

		if arguments['--cluster']:
			command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s/annotate/parse_blast.py --jmotif '%s' --nterm %s --chunk $NUM\n" % \
						( "%03d", SCRIPT_FOLDER, arguments['--jmotif'], arguments['--nterm'] )
			if arguments['--noFallBack']: command += " --noFallBack"
			pbs = open("%s/parse.sh"%prj_tree.jgene, 'w')
			pbs.write( "#!/bin/bash\n#$ -N parse-%s\n#$ -l mem=2G\n#$ -cwd\n\n%s\n" % (prj_name, command) )
			pbs.close()
			os.system( "%s -t 1-%d %s/parse.sh"%(qsub,maxFiles,prj_tree.jgene) )
		
			restart = "%s/annotate/1.3-finalize_assignments.py --reenter" % SCRIPT_FOLDER
			for opt in [ '--file', '--min1', '--min2', '--id', '--maxgaps', '--rearrangements', '--save']: 
				if arguments[opt] is not None:
					restart += " %s %s" % (opt, arguments[opt])
			for flag in ['--noclean', '--runClustering', '--runCellStatistics']:
				if arguments[flag]:
					restart += " %s" % flag

			monitor = open("%s/parse_monitor.sh"%prj_tree.jgene, 'w')
			monitor.write( "#!/bin/bash\n#$ -N monitor-%s\n#$ -l mem=2G\n#$ -cwd\n#$ -hold_jid parse-%s\n\n%s\n"%(prj_name, prj_name,restart) )
			monitor.close()
			os.system( "%s %s/parse_monitor.sh"%(qsub,prj_tree.jgene) )
			sys.exit()

		else: #do it locally

			parse_pool = Pool(arguments['--threads'])
			parse_pool.map(callParser, range(1,maxFiles+1))
			parse_pool.close()
			parse_pool.join()


	#ok, now collect all of the partial outputs and merge them
	print( "collecting information...")

	#open fasta outputs
	allV_aa	     = open ("%s/%s_allV.fa"	 % (prj_tree.aa, prj_name), "w" )
	allV_nt	     = open( "%s/%s_allV.fa"	 % (prj_tree.nt, prj_name), "w" )

	allJ_aa	     = open( "%s/%s_allJ.fa"	 % (prj_tree.aa, prj_name), "w" )
	allJ_nt	     = open( "%s/%s_allJ.fa"	 % (prj_tree.nt, prj_name), "w" )

	vj_aa	     = open( "%s/%s_goodVJ.fa"	 % (prj_tree.aa, prj_name), "w" )
	vj_nt	     = open( "%s/%s_goodVJ.fa"	 % (prj_tree.nt, prj_name), "w" )

	good_cdr3_aa = open( "%s/%s_goodCDR3.fa" % (prj_tree.aa, prj_name), "w" )
	good_cdr3_nt = open( "%s/%s_goodCDR3.fa" % (prj_tree.nt, prj_name), "w" )

	all_cdr3_aa  = open( "%s/%s_allCDR3.fa"	 % (prj_tree.aa, prj_name), "w" )
	all_cdr3_nt  = open( "%s/%s_allCDR3.fa"	 % (prj_tree.nt, prj_name), "w" )


	#also open final rearrangements tsv
	seq_stats = airr.create_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name), fields=['vj_in_frame','stop_codon','locus','c_call','junction_length','source_file','source_id','duplicate_count','length_raw','length_trimmed','indels','status','blast_identity','consensus_count','cell_id'])


	#initiate overall counters
	raw_count, total = 0, 0
	counts = {'good':0,'nonproductive':0,'indel':0,'noCDR3':0,'stop':0,'noV':0,'noJ':0,'missingNterm':0}

	dict_jcounts = Counter()
	dict_ccounts = Counter()
	dict_dcounts = Counter()
		
	c = False
	if os.path.isfile("%s/%s_C_001.txt" % (prj_tree.jgene, prj_name)):
		c = True

	d = False
	if os.path.isfile("%s/%s_D_001.txt" % (prj_tree.jgene, prj_name)):
		d = True


	#iterate over subset rearrangement files and combine
	#include generating fasta output as appropriate
	for f_ind in range(1, maxFiles+1):

		#merge partial blast hit tables
		with open( "%s/%s_jgerm_tophit.txt" % (prj_tree.tables, prj_name), "a") as table:
			with open( "%s/jtophit_%03d.txt" % (prj_tree.jgene, f_ind), "r" ) as partial:
				table.write(partial.read())

		if d:
			with open( "%s/%s_dgerm_tophit.txt" % (prj_tree.tables, prj_name), "a") as table:
				with open( "%s/dtophit_%03d.txt" % (prj_tree.jgene, f_ind), "r" ) as partial:
					table.write(partial.read())

		if c:
			with open( "%s/%s_cgerm_tophit.txt" % (prj_tree.tables, prj_name), "a") as table:
				with open( "%s/ctophit_%03d.txt" % (prj_tree.jgene, f_ind), "r" ) as partial:
					table.write(partial.read())

		#go through partial rearrangements files
		for r in airr.read_rearrangement( "%s/rearrangements_%03d.tsv"%(prj_tree.internal, f_ind) ):

			seq_stats.write( r )

			#count j/d/c gene usages
			if not r['j_call'] == "":
				dict_jcounts[ r['j_call'].split(",")[0] ] += 1
			if not r['j_call'] == "":
				dict_jcounts[ r['d_call'].split(",")[0] ] += 1
			if not r['j_call'] == "":
				dict_jcounts[ r['c_call'].split(",")[0] ] += 1

			#count statuses
			counts[ r['status'] ] += 1
			total += 1
			raw_count = int( r['sequence_id'] ) #technically, this undercounts if the last one
												# isn't in the `correct_length` interval, but I
												# don't have a better solution that isn't super
												# kludgy right now

			#ok, now do sequence output
			# start by collecting metadata for fasta def line
			def_line = ">%s" % r['sequence_id']
			if not r['v_call'] == '':          def_line += " v_call=%s"          % r['v_call']
			if not r['d_call'] == '':          def_line += " d_call=%s"          % r['d_call']
			if not r['j_call'] == '':          def_line += " j_call=%s"          % r['j_call']
			if not r['locus']  == '':          def_line += " locus=%s"           % r['locus']
			if not r['c_call'] == '':          def_line += " c_call=%s"          % r['c_call']
			if not r['status'] == '':          def_line += " status=%s"          % r['status']
#			if not r['v_identity'] == '':      def_line += " v_identity=%s"      % r['v_identity']
			if not r['junction_length'] == '': def_line += " junction_length=%s" % r['junction_length']
			if not r['junction'] == '':        def_line += " junction=%s"        % r['junction']
			if not r['junction_aa'] == '':     def_line += " junction_aa=%s"     % r['junction_aa']
			if not r['duplicate_count'] == '': def_line += " duplicate_count=%s" % r['duplicate_count']
			if not r['consensus_count'] == '': def_line += " consensus_count=%s" % r['consensus_count']
			if not r['cell_id'] == '':         def_line += " cell_id=%s"         % r['cell_id']

			#work our way up the hierarchy, putting sequences in the appropriate files
			ungapped = re.sub( "-", "", r['sequence_alignment']) #reintroduces any frameshift errors in translation
																 #  this has always been the behavior, but I wonder
																 #  if I should change/update now that I am using
																 #  proper alignments.

			if not r['status'] in ['noV', 'missingNterm']:
				allV_nt.write( "%s\n%s\n" % (def_line, ungapped) )
				allV_aa.write( "%s\n%s\n" % (def_line, Seq.Seq(ungapped).translate()) )
	
				if not r['status'] == 'noJ':
					allJ_nt.write( "%s\n%s\n" % (def_line, ungapped) )
					allJ_aa.write( "%s\n%s\n" % (def_line, Seq.Seq(ungapped).translate()) )

					if not r['status'] == 'noCDR3':
						all_cdr3_nt.write( "%s\n%s\n" % (def_line, r['junction']) )
						all_cdr3_aa.write( "%s\n%s\n" % (def_line, r['junction_aa']) )

						if r['status'] == "good":
							vj_nt.write( "%s\n%s\n" % (def_line, ungapped) )
							vj_aa.write( "%s\n%s\n" % (def_line, Seq.Seq(ungapped).translate()) )
							good_cdr3_nt.write( "%s\n%s\n" % (def_line, r['junction']) )
							good_cdr3_aa.write( "%s\n%s\n" % (def_line, r['junction_aa']) )


	#close outputs
	allV_aa.close()
	allV_nt.close()
	allJ_aa.close()
	allJ_nt.close()
	vj_aa.close()
	vj_nt.close()
	good_cdr3_aa.close()
	good_cdr3_nt.close()
	all_cdr3_aa.close()
	all_cdr3_nt.close()

	#useful number
	found = total - counts['noV'] - counts['noJ']

	#print out some statistics
	handle = open("%s/%s_jgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
	writer	= csv.writer(handle, delimiter = sep)
	keys	= sorted(dict_jcounts.keys())
	writer.writerow(["gene", "count", "percent"])
	for key in keys:
		aline = [ key, dict_jcounts[key], "%4.2f" % (dict_jcounts[key] / float(found) * 100) ]
		writer.writerow(aline)
	handle.close()

	if len(dict_ccounts) > 0:
		handle = open("%s/%s_cgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
		writer	= csv.writer(handle, delimiter = sep)
		keys	= sorted(dict_ccounts.keys())
		writer.writerow(["gene", "count", "percent"])
		for key in keys:
			aline = [ key, dict_ccounts[key], "%4.2f" % (dict_ccounts[key] / float(found) * 100) ]
			writer.writerow(aline)
		handle.close()

	if len(dict_dcounts) > 0:
		handle = open("%s/%s_dgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
		writer	= csv.writer(handle, delimiter = sep)
		keys	= sorted(dict_dcounts.keys())
		writer.writerow(["gene", "count", "percent"])
		for key in keys:
			aline = [ key, dict_dcounts[key], "%4.2f" % (dict_dcounts[key] / float(found) * 100) ]
			writer.writerow(aline)
		handle.close()

	message = "\nTotal raw reads: %d\nCorrect Length: %d\nV assigned: %d\nJ assigned: %d\nCDR3 assigned: %d\nIn-frame junction: %d\nNo indels: %d\nContinuous ORF with no stop codons: %d\n\n"  % \
								(raw_count, total, total-counts['noV'], found, found-counts['noCDR3'], found-counts['noCDR3']-counts['nonproductive'], found-counts['noCDR3']-counts['nonproductive']-counts['indel'], counts['good'])
	print( message )
	handle = open("%s/finalize_blast.log"%prj_tree.logs, "w")
	handle.write(message)
	handle.close()

	# call 1.4 if requested
	if arguments['--runClustering']:
		cmd = "%s/annotate/1.4-cluster_sequences.py" % SCRIPT_FOLDER
		for opt in [ '--file', '--min1', '--min2', '--id', '--maxgaps', '--rearrangements', '--save']: 
			if arguments[opt] is not None:
				cmd += " %s '%s'" % (opt, arguments[opt])
		if arguments['--runCellStatistics']:
			cmd += " --runCellStatistics"

		print( "Calling 1.4 with command line: %s" % cmd )
		os.system( cmd )

	#clean up!!
	oldFiles = glob.glob("%s/*txt"%prj_tree.vgene) + glob.glob("%s/*fasta"%prj_tree.vgene) +  glob.glob("%s/*txt"%prj_tree.jgene) + glob.glob("%s/*fasta"%prj_tree.jgene) + glob.glob("%s/*tsv"%prj_tree.jgene) + glob.glob("%s/lookup*"%prj_tree.internal)
	if len(oldFiles) > 0 and not arguments['--noclean']:
		[os.remove(f) for f in oldFiles]


if __name__ == '__main__':

	arguments = docopt(__doc__)
	arguments['--threads'] = int(arguments['--threads'])

	if arguments['--cluster']:
		if not clusterExists:
			sys.exit("Cannot submit jobs to non-existent cluster! Please re-run setup.sh to add support for a cluster\n")

	if arguments['--nterm'] not in ["truncate", "extend", "germline", "discard"]:
		sys.exit("--nterm must be one of ('truncate', 'extend', 'germline', 'discard')")

	#log command line
	logCmdLine(sys.argv)

	prj_tree  = ProjectFolders(os.getcwd())
	prj_name  = fullpath2last_folder(prj_tree.home)

	#load saved locus and library information
	handle = open( "%s/gene_locus.txt" % prj_tree.internal, "r")
	locus = handle.readline().strip()

	if arguments['--jmotif'] is None:
		arguments['--jmotif'] = "TGGGG"
		if "K" in locus or "L" in locus: #it's a light chain!
			if "H" in locus: #need both motifs
				arguments['--jmotif'] = "(TGGGG|TT[C|T][G|A]G)"
			else:
				arguments['--jmotif'] = "TT[C|T][G|A]G"
		elif "C" in locus and not arguments['--reenter']:
			#custom library, but default to looking for both motifs
			sys.stderr.write("Custom gene libraries used but no J motif specified; defaulting to human heavy+light...\n")
			arguments['--jmotif'] = "(TGGGG|TT[C|T][G|A]G)"

	main()
