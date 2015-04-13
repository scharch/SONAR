#!/usr/bin/python

"""
2.3-lineage-CDR3_tracebacks.py

This script compares CDR3 sequences from deep sequencing data to known CDR3 
      sequences to help identify members of a specific lineage.
                                       
Usage: 2.3-lineage-CDR3_tracebacks.py -nat natives.fa -d data.fa [-v 1-2 -j 1 -i 90 -c 90] -limits -cluster -h

Created by Zhenhai Zhang on 2011-06-30.
Modified to current form by Chaim A Schramm on 2014-01-06.
Edited and commented for publication by Chaim A Schramm on 2015-04-13.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""
import sys, os
from zap.lineage import *


def main():
	f_ind = 1

	if cluster:
		fa_file 	= "%s/C%05d.fasta" %(prj_tree.clustal_fasta, f_ind)
		print fa_file
		cmd_file 	= "%s/submitall.sh" %(prj_tree.clustal_job)
		fa_handle 	= open(fa_file, "w")
		cmd_handle 	= open(cmd_file, "w")

		for nat in natives.values():
			fa_handle.write(">%s\n" % nat.seqid)
			fa_handle.write("%s\n" % nat.seq)
	

	good = 0
	input_handle = open("%s" % dataFile, "rU")
	for ind, entry in enumerate(SeqIO.parse(input_handle, "fasta")):
		if (not germFlag) or (re.search(germline,entry.description) is not None):

			entry.seq = entry.seq[ trim[0], trim[1] ] #trim to Kabat

			if cluster:
				fa_handle.write(">%s %s\n" %(entry.id, entry.description))
				fa_handle.write("%s\n\n" %entry.seq)
			else:
				aline = [entry.id]
				for nat in natives:
					tempName = do_clustalw(nat, entry)
					percentID, percentCover = parse_pair_clustal("%s.aln"%tempFile, nat.seq_len,1) #this is one-sided coverage!
					if percentID >= ident and percentCover >= cover:
						aline.append(1)
					else:
						aline.append(0)

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

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	# get parameters from input
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)

	#check cluster parameter
	cluster = False
	flag = [x for x in ["c", "-c", "--c", "cluster", "-cluster", "--cluster"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		cluster = True

	# get parameters from input
	dict_args = processParas(sys.argv, nat="nativeFile", d="dataFile", limits="limits", v="vgene", j="jgene", i="ident", c="cover")
        default = { 'limits' : "cw", 'dataFile' : "", 'vgene' : "", 'jgene' : "", 'ident' : 90, 'cover' : 90 }
        nativeFile, dataFile, vgene, jgene, ident, cover = getParasWithDefaults(dict_args, default, "nativeFile", "dataFile", "vgene", "jgene", "ident", "cover")

	if not limits.lower() in ["cw", "imgt", "kabat"]:
		print "Format of CDR3 sequences was not recognized. Please specify 'CW', 'IMGT', or 'Kabat'.\n"
		sys.exit(0)

	trim = [9,-3]
	if limits.lower() == "kabat":
		trim = [None,None]
	elif limits.lower() == "imgt":
		trim = [6, None]

	natives = load_fastas(nativeFile)
	for nat in natives:
		nat.seq = nat.seq[ trim[0], trim[1] ]

        germFlag = True
	toMatch = ""
        if vgene == "": 
		germFlag = False
	elif jgene == "":
		toMatch = "IGHV%s"%vgene
	else:
		toMatch = "IGHV%s.*IGHJ%s"%(vgene,jgene)

	germline = re.compile(toMatch)

	prj_tree 	 = ProjectFolders(os.getcwd())
	prj_name 	 = fullpath2last_folder(prj_tree.home)
	
	if dataFile == "":
		dataFile = "%s/%s_goodCDR3.fa" % (prj_tree.nt, prj_name)

	main()

