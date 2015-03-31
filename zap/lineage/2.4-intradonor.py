#!/usr/bin/env python
# encoding: utf-8
"""
1.55-cross_donor_single_script_versionCAS.py -l is_light -npf num_per_file -g germline_gene -a use_all

Created by Zhenhai Zhang on 2011-07-12.
Modified to compress into a single script and many updates by Chaim A Schramm 2014-01-09.
Copyright (c) 2011,2014 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.55-cross_donor_single_script_versionCAS.py -l is_light -npf num_per_file -g germline_gene -a use_all

'use_all' specifies to start with the 5trim file rather than the VDJtrim

"""
import sys, os
from mytools import *
import dendropy

# zhenhai 

# mock variables, need to be masked after debugging
#HOME_FOLDER = os.getcwd()
#HOME_FOLDER = "%s/tmp" %HOME_FOLDER


# *** Zhenhai added the following block for dendropy:  2011-12-07
def contains_all_native(leafs):
	return len(leafs.intersection(native_set)) == native_len 

def retrieve_all_segregated(mytree):
	vrc01 = mytree.find_node_with_taxon_label(DICT_PHYLO_NITIVES[is_light][0])
	
	leaf_taxons = set()
	parent_node = vrc01.parent_node
	while parent_node != None:
		leafs = parent_node.leaf_nodes()
		leaf_taxons = set([str(leaf.taxon) for leaf in leafs])
		if contains_all_native(leaf_taxons):
			break
		else:
			parent_node = parent_node.parent_node
	
	#leaf_taxons = leaf_taxons.difference(native_set)
	leaf_taxons -= (native_set | set([germline_gene])) #need to remove germline for rare cases where the root node is multifurcated
	
	return leaf_taxons


def write_phylo_pbs(f):
	head = "PHY%05d" %f

	#make sure neighbor prompts us for output files
	os.system("touch outfile")
	os.system("touch outdir")

        #write a script for neighbor
	neighbor = open("%s/%s.in"%(prj_tree.clustal_fasta, head), "w")
        neighbor.write("%s/PHY%05d.mat\n"%(prj_tree.clustal_fasta,f))
        neighbor.write("F\n")
        neighbor.write("%s.out\n"%head)
        neighbor.write("2\n")
        neighbor.write("Y\n")
        neighbor.write("F\n")
        neighbor.write("%s.dnd\n"%head)
        neighbor.close()


	pbs_file = "%s/%s.sh" %(prj_tree.clustal_pbs, head)
	handle = open(pbs_file, "w")
	handle.write("#!/bin/bash\n")
	handle.write("#$ -N %s\n" %head)
	handle.write("#$ -l mem=4G,time=10:00:00\n")
	handle.write("#$ -cwd\n")
	handle.write("/ifs/home/c2b2/bh_lab/cs3037/bin/clustalo-1.1.0-linux-64 -i %s/%s.fa -o %s/%s.aln --distmat-out=%s/%s.mat --full\n" %(prj_tree.clustal_fasta, head, prj_tree.clustal_fasta, head, prj_tree.clustal_fasta, head))
	#handle.write("/ifs/home/c2b2/bh_lab/shares/clustalw/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2 -infile=%s/%s -align -tree -outputtree=nj\n" %(prj_tree.clustal_fasta, f))
	handle.write("/ifs/home/c2b2/bh_lab/cs3037/bin/neighbor < %s/%s.in"%(prj_tree.clustal_fasta,head))
	handle.close()
	return pbs_file


def get_smallest_subtree_non_natives(trees):
	for tree in trees:
		if len(tree) < native_len:
			pass;
		else:
			set_tree = set(tree)
			if len(set_tree & native_set) == native_len:
				return set_tree - native_set
	return None

def split_ids(l, num):
	llen = len(l)
	starts = range(0, llen, num)

	for start in starts:
		end = start + num
		yield set(l[start: end])




def main():
	global num_per_file, done
	# parse dnd file and get all reads clustered with native antibodies
	dnd_files, total, set_reads = sorted(glob.glob("%s/PHY*.dnd" %prj_tree.clustal_fasta)), 0, set()
	if len(dnd_files) > 0:

		#added 20120928 by CAS to handle allele wildcard characters
		os.system("perl -i -pe \"s/%s/'%s'/\" %s/PHY*.dnd" %(germline_gene.replace("*","\*"), germline_gene, prj_tree.clustal_fasta))
	
		# Zhenhai added the following block for dendropy:  2011-12-07
		for dnd_file in dnd_files:
			print "processing %s..." %dnd_file,
			tree = dendropy.Tree.get_from_path(dnd_file, "newick")
	
			outgroup_node = tree.find_node_with_taxon_label(germline_gene)
			tree.reroot_at_node(outgroup_node.parent_node)
			
			leaf_set = retrieve_all_segregated(tree)
		
			set_reads = set_reads | leaf_set
			total += len(leaf_set)
			
			if done: #kind of slow and no need to check if we've already proven it wrong
				all_taxons = set([str(leaf.taxon) for leaf in tree.leaf_nodes()])
				if len(leaf_set) < len(all_taxons) - len(native_set) - 1:
					done = False #some sequences were outside the subtree, so we haven't converged

			print "current: %d; total %d" %(len(leaf_set), total)
	
		all_reads = list(set_reads)
		if len(all_reads) == 0:
			print "NO positive sequences found --stopped!"
			sys.exit(0)

	if (len(dnd_files) == 0) or (not done):

		infasta = "%s/%s_VJtrim.fa" %(prj_tree.filtered, prj_name)
		if use_all > 0:
			infasta = "%s/%s_5trim.fa" %(prj_tree.filtered, prj_name)

		if len(dnd_files) == 0:
			print "No DND files found, starting a new XD analysis from scratch..."
			dict_reads = load_fasta_dict(infasta)
			all_reads = dict_reads.keys()
		else:
			dict_reads = load_fastas_in_set(infasta, all_reads)
	
		print "Read per file: %d!" %num_per_file
	
		# delete all previous pbs jobs in the prj/3-clustal/pbs
		os.chdir(prj_tree.clustal_pbs)
		infiles = glob.glob("*.sh")
		for infile in infiles:
			os.remove(infile)
		
		# delete all previous fasta files in the prj/3-clustal/fasta
		os.chdir(prj_tree.clustal_fasta)
		os.system("mkdir -p archives")
		oldFiles = glob.glob("archives/*")
		for old in oldFiles:
			os.remove(old)
		lastRound = glob.glob("PHY*")
		for infile in lastRound:
			os.rename(infile,"archives/%s"%infile)
	
	
		# randomize and output
		random.shuffle(all_reads)
		for f_ind, ids_set in enumerate(split_ids(all_reads, num_per_file)):
			outfile = "PHY%05d.fa" %f_ind
			handle = open(outfile, "w")
		
			# write native antibodies
			for myseq in natives:
				handle.write(">%s\n" %myseq.seq_id)
				handle.write("%s\n" %myseq.seq)
		
			for read in ids_set:
				handle.write(">%s\n" %dict_reads[read].description)
				handle.write("%s\n" %dict_reads[read].seq.tostring())
		
			handle.close()

			# write qsub script and que the job, then sleep for couple second
			pbs_file = write_phylo_pbs(f_ind)
			os.system("qsub %s" %pbs_file)

		
	elif done:
		output = open("%s/cross-donor-positives.txt"%prj_tree.data, "w")
		for read in all_reads:
			output.write("%s\n"%read)
		output.close()
		print "Tree has converged with %d reads! Output in analysis/data...\n" % len(all_reads)

	

if __name__ == '__main__':

	# get parameters from input
	if len(sys.argv) < 9 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, l="is_light", npf="num_per_file", g="germline_gene", a="use_all")
	is_light, num_per_file, germline_gene, use_all = getParas(dict_args, "is_light", "num_per_file", "germline_gene", "use_all")


	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	= fullpath2last_folder(prj_tree.home)
	nat_db		= dict_nat_db[is_light]
	
	natives		= load_list_fastas(nat_db, DICT_PHYLO_NITIVES[is_light])
	
	#added by CAS 2012/09/28
	natives.append(get_germline(germline_gene, is_light))

	native_set	= set(DICT_PHYLO_NITIVES[is_light])
	native_len  = len(DICT_PHYLO_NITIVES[is_light])

	done = True #keep track of convergence; if this is still true after processing all files, we have converged
	
	main()

