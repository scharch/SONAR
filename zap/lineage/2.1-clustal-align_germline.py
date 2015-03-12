#!/usr/bin/env python
# encoding: utf-8
"""
1.1p-split4clustal.py

Created by Zhenhai Zhang on 2011-06-30.
Copyright (c) 2011 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.2-split4clustal.py -npf num_per_file -l is_light -b use_bad


"""
import sys, os
from mytools import *

def write_a_pbs(f_ind, fa_file):
	handle = open("%s/S%07d.sh" %(prj_tree.clustal_pbs, f_ind), "w")
	handle.write("#!/bin/bash\n")
	handle.write("#$ -N S%07d\n" %(f_ind))
	handle.write("#$ -l mem=1G,time=30:00:00\n")
	handle.write("#$ -cwd\n")
	
	# Usage: 1.2-clustalw2.py -l is_light -i infile -o outfile -d data_file
	handle.write("clustalw2.py -l %d -i %s -o %s/S%07d.txt"  %(is_light, fa_file, prj_tree.clustal_data, f_ind))
	handle.close()
	

def main():
	f_ind = 1
	fa_file 	= "%s/S%07d.fasta" %(prj_tree.clustal_fasta, f_ind)
	print fa_file
	cmd_file 	= "%s/submitall.sh" %(prj_tree.clustal_job)
	fa_handle 	= open(fa_file, "w")
	cmd_handle 	= open(cmd_file, "w")

	if use_bad:
		input_handle = open("%s/%s_allJ.fa" %(prj_tree.filtered, prj_name), "rU")
	else:
		input_handle = open("%s/%s_VJtrim.fa" %(prj_tree.filtered, prj_name), "rU")

	for ind, entry in enumerate(SeqIO.parse(input_handle, "fasta")):
		fa_handle.write(">%s %s\n" %(entry.id, entry.description))
		fa_handle.write("%s\n\n" %entry.seq.tostring())
		
		if (ind + 1) % num_per_file == 0:
			write_a_pbs(f_ind, fa_file)
			cmd_handle.write("qsub %s/S%07d.sh\n" %(prj_tree.clustal_pbs, f_ind))
			#cmd_handle.write("sleep 5\n\n")
			f_ind += 1
			fa_handle.close()
			fa_file 	= "%s/S%07d.fasta" %(prj_tree.clustal_fasta, f_ind)
			print fa_file
			fa_handle = open(fa_file, "w")
	
	# last file 	
	if (ind + 1) % num_per_file > 0:
		write_a_pbs(f_ind, fa_file)
		cmd_handle.write("qsub %s/S%07d.sh\n" %(prj_tree.clustal_pbs, f_ind))
	
	# submit all jobs
	cmd_handle.close()
	os.chdir(prj_tree.clustal_job)
	os.system("./submitall.sh > submit_all.txt &")

if __name__ == '__main__':

	# get parameters from input
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)

	dict_args 	 = processParas(sys.argv, npf="npf", l="is_light", b="use_bad")
	num_per_file, is_light, use_bad = getParas(dict_args, "npf", "is_light", "use_bad")
	
	prj_tree 	 = ProjectFolders(os.getcwd())
	prj_name 	 = fullpath2last_folder(prj_tree.home)
	
	main()

