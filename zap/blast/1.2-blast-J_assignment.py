#!/usr/bin/python

"""
1.2-blast-J_assignment.py

This script parses the BLAST output from 1.1-blast-V_assignment.py. For reads
      for which a V assignment was successfully made, the section of the read
      3' to the V gene is extracted and sent to the cluster for BLAST assignment
      of the J gene.

Usage: 1.2-blast-J_assignment.py -locus <0|1|2|3|4> -lib path/to/library.fa
                                 -dlib path/to/dlibrary.fa -h

    All options are optional, see below for defaults.
    Invoke with -h or --help to print this documentation.

    locus	0: heavy chain / 1: kappa chain / 2: lambda chain / 3: kappa OR
                   lambda / 4: custom library (supply -lib and optionally -dlib)
                   Default = 0

Created by Zhenhai Zhang on 2011-04-14.
Modified by Chaim A Schramm 2013-07-03 to include j assignment.
Modified by Chaim A Schramm 2014-03-25 to not swamp RAM when processing Illumina
    data.
Edited and commented for publication by Chaim A Schramm on 2015-02-09.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""

import sys
import os
from zap.blast import *


def get_top_hit(resultFile, writer):

	dict_germ_aln=dict()
	
	for my_alignment, row, others, second_match in generate_blast_top_hists(resultFile):
		qid, sid 	= 	my_alignment.qid, my_alignment.sid
		aline 		= 	row + [my_alignment.strand]
		
		dict_germ_aln[qid] = my_alignment

		if len(others)>0:
			aline.append(",".join(others))
		writer.writerow(aline)

		if len(second_match)>0:
			writer.writerow(second_match + [my_alignment.strand])
			#if shorter match is toward 3', change alignment end to get J properly
			if my_alignment.strand=="+" and int(second_match[7]) > my_alignment.qend:
				my_alignment.qend = int(second_match[7])
			if my_alignment.strand=="-" and int(second_match[6]) < my_alignment.qstart:
				my_alignment.qstart = int(second_match[6])
		
		if my_alignment.sid not in dict_germ_count:
			dict_germ_count[my_alignment.sid] = 0
		dict_germ_count[my_alignment.sid] += 1
		
	return dict_germ_aln



def main():
	
	print "curating 5'end and strand...."

	# cut nucleotide sequences from 5'end alignment to germline
	total, good, f_ind = 0, 0, 1

	writer = csv.writer(open("%s/%s_vgerm_tophit.txt" %(prj_tree.tables, prj_name), "w"), delimiter = sep)
	writer.writerow(PARSED_BLAST_HEADER)
	
	while os.path.isfile("%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind)):

		split_fasta = "%s/%s_%03d.fasta" %(prj_tree.jgene, prj_name, f_ind)
		fasta_handle = open(split_fasta, "w")

		# parse blast output
		dict_germ_aln = get_top_hit("%s/%s_%03d.txt" % (prj_tree.vgene, prj_name, f_ind), writer)
	
		# process each sequence
		for entry in SeqIO.parse("%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind), "fasta"):

			try:
				#to match BLAST, which also kills zero-padding
				seq_id = str(int(entry.id))
			except:
				seq_id = entry.id

			total += 1
			if seq_id in dict_germ_aln:
				if dict_germ_aln[seq_id].strand == "+":
					entry.seq = entry.seq[ dict_germ_aln[seq_id].qend : ]
				else:
					entry.seq = entry.seq[ : dict_germ_aln[seq_id].qstart]
					entry = entry.reverse_complement()

				if len(entry.seq) > 30: #can probably be 50...
					fasta_handle.write(">%s\n%s\n" % (entry.id,entry.seq))
					good += 1

		fasta_handle.close()
		f_ind += 1
		
		print "%d done, %d good..." %(total, good)

	f_ind -= 1 #had to go 1 extra to break while loop, now reset to actual number of files

	#print statistics
	handle = open("%s/%s_vgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
	writer 	= csv.writer(handle, delimiter = sep)
	keys 	= sorted(dict_germ_count.keys())
	
	writer.writerow(["gene", "count", "percent"])
	for key in keys:
		aline = [ key, dict_germ_count[key], "%4.2f" % (dict_germ_count[key] / float(good) * 100) ]
		writer.writerow(aline)
		
	
	#print log message
	handle = open("%s/1.2.log" % prj_tree.logs, "w")
	handle.write("total: %d; good: %d\n" %(total, good))
	handle.close()
	
	# write pbs files and auto submit shell script
	command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLASTALL % (BLAST_J_OPTIONS, library, 
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_$NUM.txt"   % (prj_tree.jgene, prj_name)) )
	pbs = open("%s/jblast.sh"%prj_tree.jgene, 'w')
	pbs.write( PBS_STRING%(f_ind, "%s-jBlast"%prj_name, "500M", "10:00:00", "%s 2> %s/%s_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
	pbs.close()
	os.system("qsub %s/jblast.sh"%prj_tree.jgene)


	if os.path.isfile(dlib):
		command = "NUM=`printf \"%s\" $SGE_TASK_ID`\n%s" % ( "%03d", CMD_BLASTALL % (BLAST_J_OPTIONS, dlib, 
									      "%s/%s_$NUM.fasta" % (prj_tree.jgene, prj_name),
									      "%s/%s_D_$NUM.txt"   % (prj_tree.jgene, prj_name)) )
		pbs = open("%s/dblast.sh"%prj_tree.jgene, 'w')
		pbs.write( PBS_STRING%(f_ind, "%s-dBlast"%prj_name, "500M", "10:00:00", "%s 2> %s/%s_D_$NUM.err"%(command, prj_tree.jgene, prj_name)) )
		pbs.close()
		os.system("qsub %s/dblast.sh"%prj_tree.jgene)


if __name__ == '__main__':
	
	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	#get parameters from input
	dict_args = processParas(sys.argv, locus="locus", library="library", dlib="dlib")
	locus, library, dlib = getParasWithDefaults(dict_args, dict(locus=0, library="", dlib=""), "locus", "library", "dlib")

	#load library
	if locus < 4:
		library = dict_jgerm_db[locus]
		if locus == 0:
			dlib = DH_DB
	elif os.path.isfile(library):
		pass
	else:
		print "Can't find custom library file!"
		sys.exit(1)

	prj_tree        = ProjectFolders(os.getcwd())
	prj_name        = fullpath2last_folder(prj_tree.home)

	dict_germ_count	= dict()
	
	
	main()

