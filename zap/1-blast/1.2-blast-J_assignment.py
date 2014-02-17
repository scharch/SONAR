#!/usr/bin/env python
# encoding: utf-8
"""
1.1-parse_blast_result.py

Created by Zhenhai Zhang on 2011-04-14.
Modified by CA Schramm 2013-07-03 to include j assignment.
Copyright (c) 2011, 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.1-parse_blast_result.py -l [0: heavy; 1: kappa; 2: lambda; 3: k/l] -npf 50000


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
		
	

def get_top_hit(folder):

	#qstart_list 			= [0] * (MAX_GOOD_LEN + 1)
	#qend_list				= [0] * (MAX_GOOD_LEN + 1)
	#sstart_list				= [0] * (MAX_GOOD_LEN + 1)
	#send_list				= [0] * (MAX_GOOD_LEN + 1)
	
	dict_strand_count 		= {"+" : 0.0, "-" : 0.0}
	dict_germ_count			= dict()
	
	db_name = fullpath2last_folder(folder)
	writer, total = csv.writer(open("%s/%s_vgerm_tophit.txt" %(prj_tree.data, prj_name), "w"), delimiter = sep), 0
		
	writer.writerow(PARSED_BLAST_HEADER)
	
	
	for my_alignment, row, others in generate_blast_top_hists(folder, prj_name):
		qid, sid 	= 	my_alignment.qid, my_alignment.sid
		aline 		= 	row + [my_alignment.strand]
		
		if db_name == "germ":
			aline.append(dict_germ[sid].seq_len)
			# record germline alignment
			dict_germ_aln[qid] = my_alignment
		elif db_name == "native":
			sid = sid.upper()
			aline.append(dict_native[sid].seq_len)
		else:
			print "DATABASE NAME ERROR!"
			sys.exit(0)		

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
		
		total		+= 	1
		

	#write_list2file(qstart_list, "%s_qstart" %db_name, total)
	#write_list2file(sstart_list, "%s_sstart"  %db_name, total)
	#write_list2file(qend_list, "%s_qend" %db_name, total)
	#write_list2file(send_list, "%s_send" %db_name, total)
	
	write_dict2file(dict_strand_count, total, ["strand", "count", "percent"])
	write_dict2file(dict_germ_count, total, ["subject", "count", "percent"])	


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
		
	# parse germline alignment
	get_top_hit(prj_tree.germ)
	
	#outfile = "%s/%s_5trim.fa" %(prj_tree.filtered, prj_name)
	#handle = open(outfile, "w")
	
	
	print "curating 5'end and strand...."
	# cut nucleotide sequences from 5'end alignment to germline
	total, good, f_ind = 0, 0, 1
	split_fasta = "%s/vcut_%06d.fasta" %(prj_tree.split, f_ind)
	fasta_handle = open(split_fasta, "w")

	for entry in SeqIO.parse(open("%s/%s.fasta" %(prj_tree.filtered, prj_name), "rU"), "fasta"):
		seq_id = str(int(entry.id))
		if seq_id in dict_germ_aln:
			entry.description = dict_germ_aln[seq_id].sid
			entry.seq = cut_five_end_blast(entry.seq, dict_germ_aln[seq_id]).tostring()
			align_len = abs(dict_germ_aln[seq_id].qend - dict_germ_aln[seq_id].qstart) + 1
			entry.seq = entry.seq[ align_len : ]

			if len(entry.seq) > 30: #can probably be 50...
				fasta_handle.write(">%s\n%s\n" % (entry.id,entry.seq))
				good += 1

			if good % num_per_file == 0 and good>0:
				fasta_handle.close()
				submit_pbs(split_fasta, f_ind)
				f_ind += 1
				split_fasta 	= "%s/vcut_%06d.fasta" 	%(prj_tree.split, f_ind)
				fasta_handle	= open(split_fasta, "w")	
			
			'''
			if is_light == 0:  # human heavy chain
				has_wgxg, wgxg_start, wgxg_end = has_pat(entry.seq, pat_nuc_wgxg)
				has_ast, ast_start, ast_end = has_pat(entry.seq, pat_nuc_ast)
				has_vss, vss_start, vss_end = has_pat(entry.seq, pat_nuc_vss)
				if has_ast and has_wgxg and (ast_start > wgxg_end):
					entry.seq = entry.seq[ : ast_start]
				
				elif not has_ast and has_vss and has_wgxg and (vss_start > wgxg_end):
					entry.seq = entry.seq[ : vss_end]
				
			if len(entry.seq) > 50: #kludge for Xenomouse light chain; change back to: MIN_ANTIBODY_LEN:
				handle.write(">%s %s\n" %(entry.id, entry.description))
				handle.write("%s\n" %entry.seq)
				good += 1
			'''

		total += 1
		
		if total % (10 ** 5) == 0:
			print "%d done, %d good..." %(total, good)

	fasta_handle.close()
	submit_pbs(split_fasta, f_ind)	
	print "%d done, %d good..." %(total, good)


if __name__ == '__main__':
	
	# get parameters from input
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, l="is_light", npf="num_per_file")
	is_light, num_per_file = getParas(dict_args, "is_light", "num_per_file")	

	prj_tree       =  ProjectFolders(os.getcwd())
	prj_name       =  fullpath2last_folder(prj_tree.home)

	germ_db, j_db  =  dict_germ_db[is_light], dict_jgerm_db[is_light]
	dict_germ      =  load_fastas(germ_db)
	dict_j         =  load_fastas(j_db)
	dict_germ_aln  =  dict()
	
	main()

