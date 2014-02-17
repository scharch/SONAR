#!/usr/bin/env python
# encoding: utf-8
"""
1.0-preliminary.py

Created by Zhenhai Zhang on 2011-04-12.
Copyright (c) 2011 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.0-preliminary.py -nl min_len -ml max_len -npf reads_per_file -k [0: Heavy; 1: Light] -hq [1:has quals; 0: noquals]
min_len		Minimum length for read filtering (inclusive)
max_len		Maximum length for read filtering (inclusive)
k = 0: heavy chain
k = 1: light chain

"""
# --BEGIN-- general package
import sys
import os

from mytools import *
# -- END -- general package

# global variables
read_len_count, total, total_good, f_ind = [0] * MAX_READ_LEN, 0, 0, 1


def main():
	global read_len_count, total, total_good, f_ind
		
	
	# file name
	merged_file 	= "%s/%s.fasta" 			%(folder_tree.filtered, prj_name)
	split_fasta 	= "%s/%s_%06d.fasta" 		%(folder_tree.split, prj_name, f_ind)
	split_qual 		= "%s/%s_%06d.qual" 		%(folder_tree.split, prj_name, f_ind)
	read_id_map		= "%s/%s_454_rid.txt"		%(folder_tree.data, prj_name)
	
	merge_handle, fasta_handle, id_writer = open(merged_file, "w"), open(split_fasta, "w"), csv.writer(open(read_id_map, "w"), delimiter = sep)
	if has_qual:
		qual_handle 	=  open(split_qual, "w")
		qual_generater 	= generate_quals_folder(folder_tree.original)
		
	id_writer.writerow(["454_id", "read_id"])
	
	
	for myseq in generate_read_fasta_folder(folder_tree.original):
		total += 1
		
		read_len_count[myseq.seq_len] += 1
		
		if min_len <= myseq.seq_len <= max_len:
			total_good += 1
			
			# record off machine id to read id
			id_writer.writerow([myseq.seq_id, total_good])
			
			
			
			merge_handle.write(">%08d\n" 	%total_good)
			merge_handle.write("%s\n" 	%myseq.seq)
			
			fasta_handle.write(">%08d\n"  %total_good)
			fasta_handle.write("%s\n" 	%myseq.seq)
			
			if has_qual:
				myqual = qual_generater.next()
				qual_handle.write(">%08d\n" 	%total_good)
				qual_handle.write("%s\n" 		%" ".join(map(str, myqual.qual_list)))
			
			
			if total_good % read_per_file == 0:
				
				fasta_handle.close()
				
				
				f_ind += 1
				split_fasta 	= "%s/%s_%06d.fasta" 	%(folder_tree.split, prj_name, f_ind)
				fasta_handle	= open(split_fasta, "w")
				
				print parse_name(split_fasta)
				
				if has_qual:
					qual_handle.close()				
					qual_handle = open(split_qual, "w")
					split_qual 	= "%s/%s_%06d.qual" 	%(folder_tree.split, prj_name, f_ind)
		
		if total % 10 ** 4 == 0:
			print "%d processed, %d good" %(total, total_good)
	
	print "%d processed, %d good" %(total, total_good)
	handle = open("%s/1-split.log" %folder_tree.logs, "w")
	handle.write("total: %d; good: %d; percentile: %f\n" %(total, total_good, float(total_good)/total * 100))
	
	# write pbs files and auto submit shell script
	write_pbs_file(folder_tree, is_light)
	write_auto_submit_file(folder_tree)
	
	# write length file
	#write_read_length(folder_tree, read_len_count, total)
	
	# TEMPORARILY disable the auto submission 
	os.chdir(folder_tree.jobs)
	os.system("./submit_all.sh > submit_all.txt &")
	
	

if __name__ == '__main__':

	# get parameters from input
	if len(sys.argv) < 9 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, nl="min_len", ml="max_len", npf="read_per_file", k="is_light", hq="hasqual")
	min_len, max_len, read_per_file, is_light, has_qual = getParas(dict_args, "min_len", "max_len", "read_per_file", "is_light", "hasqual")
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	folder_tree = ProjectFolders(prj_folder)
	folder_tree = create_folders( prj_folder )
	
	
	prj_name = prj_folder[prj_folder.rindex("/") + 1 :]

	main()

