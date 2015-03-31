#!/usr/bin/env python

"""

                                       
Created by Chaim A Schramm on 2014-01-06. Modified from 1.2-split4clustal.py
Copyright (c) 2011,2014 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.71-CDR3_tracebacks.py -npf num_per_file -nat natives.fa -d data.fa [-v 1-2 -j 1 -i 90 -c 90]


"""
import sys, os
from mytools import *

def write_a_pbs(f_ind, fa_file):
	handle = open("%s/C%05d.sh" %(prj_tree.clustal_pbs, f_ind), "w")
	handle.write("#!/bin/bash\n")
	handle.write("#$ -N C%05d\n" %(f_ind))
	handle.write("#$ -l mem=1G,time=10:00:00\n")
	handle.write("#$ -cwd\n")
	
	# Usage: clustalw2_CAS.py -n numNat -i infile -o outfile
	handle.write("clustalw2_CAS.py -n %d -i %s -o %s/C%05d.txt -id %d -cov %d\n"  %(len(natives), fa_file, prj_tree.clustal_data, f_ind, ident, cover))
	handle.close()
	

def main():
	f_ind = 1
	fa_file 	= "%s/C%05d.fasta" %(prj_tree.clustal_fasta, f_ind)
	print fa_file
	cmd_file 	= "%s/submitall.sh" %(prj_tree.clustal_job)
	fa_handle 	= open(fa_file, "w")
	cmd_handle 	= open(cmd_file, "w")

	for nat in natives.values():
		fa_handle.write(">%s\n" % nat.seq_id)
		fa_handle.write("%s\n" % nat.seq[9:-3]) #trim CDRH3 to Kabat defintion (needs to be changed for light chains)
	
	good = 0
	input_handle = open("%s" % dataFile, "rU")
	for ind, entry in enumerate(SeqIO.parse(input_handle, "fasta")):
		if (not germFlag) or (re.search(germline,entry.description) is not None):
			fa_handle.write(">%s %s\n" %(entry.id, entry.description))
			fa_handle.write("%s\n\n" %entry.seq.tostring()[9:-3]) #trim to Kabat
			good +=1

		if good == num_per_file:
			good = 0
			write_a_pbs(f_ind, fa_file)
			cmd_handle.write("qsub %s/C%05d.sh\n" %(prj_tree.clustal_pbs, f_ind))
			#cmd_handle.write("sleep 5\n\n")
			f_ind += 1
			fa_handle.close()
			fa_file 	= "%s/C%05d.fasta" %(prj_tree.clustal_fasta, f_ind)
			print fa_file
			fa_handle = open(fa_file, "w")
			for nat in natives.values():
				fa_handle.write(">%s\n" % nat.seq_id)
				fa_handle.write("%s\n" % nat.seq[9:-3])
	
	# last file 	
	if good > 0:
		write_a_pbs(f_ind, fa_file)
		cmd_handle.write("qsub %s/C%05d.sh\n" %(prj_tree.clustal_pbs, f_ind))
	
	# submit all jobs
	cmd_handle.close()
	os.chdir(prj_tree.clustal_job)
	os.system("./submitall.sh > submit_all.txt &")

if __name__ == '__main__':

	# get parameters from input
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)

	dict_args 	 = processParas(sys.argv, npf="npf", nat="nativeFile", d="dataFile", v="vgene", j="jgene", i="ident", c="cover")
        default = { 'vgene' : "", 'jgene' : "", 'ident' : 90, 'cover' : 90 }
        num_per_file, nativeFile, dataFile, vgene, jgene, ident, cover = getParasWithDefaults(dict_args, default, "npf", "nativeFile", "dataFile", "vgene", "jgene", "ident", "cover")

        germFlag = True
        if vgene == "": germFlag = False

	print "Thresholds are %d%% identity and %d%% (one-sided) coverage." %(ident,cover)
	print "Germline regex is: \"IGHV%s.*IGHJ%s\""%(vgene,jgene)
	germline = re.compile("IGHV%s.*IGHJ%s"%(vgene,jgene))

	natives = load_fastas(nativeFile)
	
	prj_tree 	 = ProjectFolders(os.getcwd())
	prj_name 	 = fullpath2last_folder(prj_tree.home)
	
	main()

