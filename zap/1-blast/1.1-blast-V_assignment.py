#!/usr/bin/python

"""
1.1-blast-V_assignment.py

This script looks for raw NGS data files in the "0-original" folder and parses
      them into manageable chunks for searching with BLAST on the cluster. Each
      raw read that meets the length restrictions is assigned a unique 8 digit 
      serial number and submitted to the cluster for BLASTing against the 
      desired V gene germline library.
      Distinctions between raw data files will not be maintained, although they
      can be reconstructed using the id table. The name of the current folder
      is assigned as the project name, which will be used to identify all
      output files created. Resource requests for the cluster are calibrated to
      groups of 50K sequences, so that number is hard-coded below.

Usage: 1.1-blast-V_assignment.py -minl min_len -maxl max_len -locus <0|1|2|3|4>
                                 -qual <0|1|2> -lib path/to/library.fa -h -f

    All options are optional, see below for defaults.
    Invoke with -h or --help to print this documentation.

    minl	Minimum length for read filtering (inclusive). Default = 300
    maxl	Maximum length for read filtering (inclusive). Default = 600.
    locus	0: heavy chain / 1: kappa chain / 2: lambda chain / 3: kappa OR
                   lambda / 4: custom library (supply -lib)
                   Default = 0
    qual 	0: noquals/fasta only / 1: has quals (454) / 2: fastq (Illumina).
                   Default = 0 (will fail if reads are in FastQ format)
		   NOTE: quals are not currently processed, but this flag is
		         still required for correctly reading FastQ files!
    lib  	location of file containing custom library (e.g. for use with
                   non-human genes)
    f 	 	forcing flag to overwrite existing working directories.

Created by Zhenhai Zhang on 2011-04-12.
Edited and commented for publication by Chaim A Schramm on 2014-12-22.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""


import sys
import os
from zap.blast import *


# global variables
total, total_good, f_ind = 0, 0, 1


def main():

	global total, total_good, f_ind
		
	# open initial output files
	fasta 	=            open("%s/%s_%06d.fasta"  % (folder_tree.vgene,  prj_name, f_ind), 'w')
	id_map	= csv.writer(open("%s/%s_454_rid.txt" % (folder_tree.tables, prj_name),        'w'), delimiter=sep)

	id_map.writerow(["454_id", "read_id"])
	

	#if we decide to use quals for something, can add this block back in
	'''
	if has_qual > 0:
		qual =       open("%s/%s_%06d.qual"   % (folder_tree.vgene,    prj_name, f_ind), 'w')
	if has_qual == 1:
		qual_generater 	= generate_quals_folder(folder_tree.home)
	'''


	#iterate through sequences in all raw data files
	for myseq, myqual in generate_read_fasta_folder(folder_tree.home, has_qual):
		total += 1
			
		if min_len <= myseq.seq_len <= max_len:
			total_good += 1
			id_map.writerow([myseq.seq_id, total_good])
			fasta.write(">%08d\n%s\n" % (total_good, myseq.seq))

			#uncomment to re-implement quals
			'''
			if has_qual == 1:
				myqual = qual_generater.next()
			if has_qual > 0:
				qual.write(">%08d\n%s\n" % (total_good, " ".join(map(str, myqual.qual_list))))
			'''

			if total_good % 50000 == 0: 
				#close old output files, open new ones, and print progress message
				fasta.close()
				f_ind += 1
				fasta = open("%s/%s_%06d.fasta" % (folder_tree.vgene, prj_name, f_ind), 'w')
				print "%d processed, %d good; starting file %s_%06d" %(total, total_good, prj_name, f_ind)

				'''
				if has_qual>0:
					qual.close()
					qual = open("%s/%s_%06d.qual"%(folder_tree.vgene, prj_name, f_ind), 'w')
				'''
				
	print "TOTAL: %d processed, %d good" %(total, total_good)
	
	fasta.close()
	'''
	if has_qual>0:
		qual.close()
	'''

	#print log message
	handle = open("%s/1-split.log" % folder_tree.logs, "w")
	handle.write("total: %d; good: %d; percentile: %f\n" %(total, total_good, float(total_good)/total * 100))
	handle.close()
	
	# write pbs files and auto submit shell script
	write_pbs_file(folder_tree.vgene, locus, library)
	write_auto_submit_file(folder_tree.vgene, submit=True)
	


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

	# get parameters from input
	dict_args = processParas(sys.argv, minl="min_len", maxl="max_len", locus="locus", qual="has_qual", lib="library")
	min_len, max_len, locus, has_qual, library = getParasWithDefaults(dict_args, dict(min_len=300, max_len=600, has_qual=0, locus=0, library=""), "min_len", "max_len", "locus", "has_qual", "library")

	# create 1st and 2nd subfolders
	prj_folder  = os.getcwd()
	folder_tree = create_folders( prj_folder, force=force )
	prj_name    = prj_folder[prj_folder.rindex("/") + 1 :]

	
	main()

