#!/usr/bin/env python
# encoding: utf-8
"""
1.1-parse_blast_result.py

Created by Zhenhai Zhang on 2011-04-14.
Modified by CA Schramm 2013-07-03 to include j assignment.
Modified by Chaim Schramm 2014-03-25 to not swamp RAM when processing Illumina data.
Copyright (c) 2011, 2013, 2014 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.1-parse_blast_result.py -l [0: heavy; 1: kappa; 2: lambda; 3: k/l] 


"""
# --BEGIN-- general package
import sys
import os

from mytools import *
# -- END -- general package

def write_list2file(l, sname, total):
	writer = csv.writer(open("%s/%s_%s_stat.txt" %(prj_tree.data, prj_name, sname), "w"), delimiter = sep)
	percentiles = [(x / total) * 100 for x in l]
	#percentiles = l / float(total) * 100
	
	writer.writerow(["index", "count", "percentile"])
	for ind, count in enumerate(l):
		aline = [ind, count, percentiles[ind]]
		writer.writerow(aline)
		

def write_dict2file(d, total, header):
	filename = "%s/%s_vgerm_stat.txt" %(prj_tree.data, prj_name)
	if "strand" in header:
		filename = "%s/%s_germ_strand_stat.txt" %(prj_tree.data, prj_name)
	writer 	= csv.writer(open(filename, "w"), delimiter = sep)
	keys 	= sorted(d.keys())
	
	writer.writerow(header)
	
	for key in keys:
		aline = [key, d[key], d[key] / float(total) * 100]
		writer.writerow(aline)
		
	

def get_top_hit(resultFile, writer):

	#qstart_list 			= [0] * (MAX_GOOD_LEN + 1)
	#qend_list				= [0] * (MAX_GOOD_LEN + 1)
	#sstart_list				= [0] * (MAX_GOOD_LEN + 1)
	#send_list				= [0] * (MAX_GOOD_LEN + 1)
	
	#db_name = fullpath2last_folder(folder)
	dict_germ_aln=dict()
	
	for my_alignment, row, others in generate_blast_top_hists(resultFile):
		qid, sid 	= 	my_alignment.qid, my_alignment.sid
		aline 		= 	row + [my_alignment.strand]
		
		#if db_name == "germ":
		aline.append(dict_germ[sid].seq_len)
			# record germline alignment
		dict_germ_aln[qid] = my_alignment
		#elif db_name == "native":
		#	sid = sid.upper()
		#	aline.append(dict_native[sid].seq_len)
		#else:
		#	print "DATABASE NAME ERROR!"
		#	sys.exit(0)		

		if len(others)>0:
			aline.append(",".join(others))
		writer.writerow(aline)
		
		
		if my_alignment.sid not in dict_germ_count:
			dict_germ_count[my_alignment.sid] 	= 0
		dict_germ_count[my_alignment.sid] 		+= 1
		
		
		sstart, send = my_alignment.sstart, my_alignment.send
		if my_alignment.strand == "-":
			sstart, send = send, sstart

		#qstart_list[my_alignment.qstart]	+= 1
		#qend_list[my_alignment.qend]		+= 1		
		#sstart_list[sstart] += 1
		#send_list[send]		+= 1
		dict_strand_count[my_alignment.strand] += 1
		
		#total		+= 	1
		

	#write_list2file(qstart_list, "%s_qstart" %db_name, total)
	#write_list2file(sstart_list, "%s_sstart"  %db_name, total)
	#write_list2file(qend_list, "%s_qend" %db_name, total)
	#write_list2file(send_list, "%s_send" %db_name, total)
	
	return dict_germ_aln


def submit_pbs(fasta, index):
	pbsfile    = "%s/vc%06d.sh" %(prj_tree.pbs, index)
	jobname    = "vc%06d" % index
	outfile    = "%s/vc%06d.txt" %(prj_tree.germ, index)
	pbs_string = PBS_STRING %(jobname, outfile, BLAST_GERM_J_OPTIONS, j_db, fasta)
	
	pbshandle  = open( pbsfile, "w")
	pbshandle.write(pbs_string)
	pbshandle.close()
	
	os.system("qsub %s" %pbsfile)


def main():
		
	#outfile = "%s/%s_5trim.fa" %(prj_tree.filtered, prj_name)
	#handle = open(outfile, "w")
	
	
	print "curating 5'end and strand...."
	#os.system("rm %s/%s_vgerm_tophit.txt" %(prj_tree.data, prj_name))

	# cut nucleotide sequences from 5'end alignment to germline
	total, good, f_ind = 0, 0, 1

	writer = csv.writer(open("%s/%s_vgerm_tophit.txt" %(prj_tree.data, prj_name), "w"), delimiter = sep)
	writer.writerow(PARSED_BLAST_HEADER)
	
	while os.path.isfile("%s/%s_%06d.fasta" % (prj_tree.split, prj_name, f_ind)):

		split_fasta = "%s/vcut_%06d.fasta" %(prj_tree.split, f_ind)
		fasta_handle = open(split_fasta, "w")

		# parse germline alignment
		dict_germ_aln = get_top_hit("%s/%s_%06d.txt" % (prj_tree.germ, prj_name, f_ind), writer)
	
		for entry in SeqIO.parse(open("%s/%s_%06d.fasta" % (prj_tree.split, prj_name, f_ind), "rU"), "fasta"):

			seq_id = str(int(entry.id))
			total += 1
			if seq_id in dict_germ_aln:
				entry.description = dict_germ_aln[seq_id].sid
				entry.seq = cut_five_end_blast(entry.seq, dict_germ_aln[seq_id]).tostring()
				align_len = abs(dict_germ_aln[seq_id].qend - dict_germ_aln[seq_id].qstart) + 1
				entry.seq = entry.seq[ align_len : ]

				if len(entry.seq) > 30: #can probably be 50...
					fasta_handle.write(">%s\n%s\n" % (entry.id,entry.seq))
					good += 1

		fasta_handle.close()
		submit_pbs(split_fasta, f_ind)
		f_ind += 1
		
		print "%d done, %d good..." %(total, good)

	write_dict2file(dict_strand_count, total, ["strand", "count", "percent"])
	write_dict2file(dict_germ_count, total, ["subject", "count", "percent"])	


if __name__ == '__main__':
	
	# get parameters from input
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, l="is_light")#, npf="num_per_file")
	is_light = getParas(dict_args, "is_light")#, "num_per_file")	

	prj_tree       =  ProjectFolders(os.getcwd())
	prj_name       =  fullpath2last_folder(prj_tree.home)

	germ_db, j_db  =  dict_germ_db[is_light], dict_jgerm_db[is_light]
	dict_germ      =  load_fastas(germ_db)
	dict_j         =  load_fastas(j_db)
	#dict_germ_aln  =  dict()
	dict_strand_count 		= {"+" : 0.0, "-" : 0.0}
	dict_germ_count			= dict()
	
	
	main()

