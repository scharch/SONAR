#!/usr/bin/env python


"""
Core functions for the "ZAP" package,

Created by Zhenhai Zhang on 2011-04-05 as mytools.py
Edited and commented for publication by Chaim A Schramm on 2015-02-10.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.
"""

from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from math import log

import glob
import re
import subprocess


from numpy import mean, array, zeros, ones, nan, std, isnan

from commonVars import *
from numpy_related import *



#
# -- BEGIN -- class definition
#

class MyNode:
	def __init__(self, parent):
		self.parent 	= parent
		self.children 	= []
	def add_child(self, child):
		self.children.append(child)


class MySeq:
	def __init__(self, seq_id, seq):
		self.seq_id 	= seq_id					# sequence ID
		try:
			self.seq 	= str(seq)			# sequence in string format
		except:
			self.seq	= seq	
		self.seq_len 	= len(self.seq)				# sequence length
		
		
	
	def set_desc(self, desc):
		self.desc = desc
		
class MyQual:
	def __init__(self, qual_id, qual_list):
		self.qual_id 	= qual_id
		self.qual_list 	= qual_list
		self.qual_len	= len(qual_list)
		
class MyAlignment:
	def __init__(self, row):
		self.qid	= row[0].strip()		# query id
		self.sid	= row[1].upper().strip()	# subject id
		self.identity 	= float(row[2])			# % identity
		self.alignment 	= int(row[3])			# alignment length
		self.mismatches = int(row[4]) 			# mismatches
		self.gaps	= int(row[5])			# gap openings
		self.qstart 	= int(row[6])			# query start
		self.qend	= int(row[7])			# query end
		self.sstart 	= int(row[8])			# subject start
		self.send	= int(row[9]) 			# subject end
		self.evalue 	= float(row[10])		# e-value
		self.score	= float(row[11])		# bit score
		self.strand	= str(row[12])			# strand
		
		self.qlen	= 0
		self.slen	= 0
		
		self.real_id	= 0.0				# recaluclated identity
		self.divergence	= 0.0				# recalculated diversity
		
	def set_strand(self, s):
		self.strand = s					# setting strand
		
	def set_real_identity(self, identity):
		self.real_id = identity
		
	def set_diversity(self, divergence):
		self.divergence = divergence
		
class MyAlignmentVerbose:
	def __init__(self, row):
		self.query_id		= row[0]
		self.sbjct_id		= row[1].upper()
		self.strand		= row[2]
		self.evalue		= float(row[3])
		self.score		= float(row[4])
		self.identities		= int(row[5])
		self.gaps		= int(row[6])
		self.aln_len		= int(row[7])
		self.query_start	= int(row[8])
		self.query_end		= int(row[9])
		self.query_len		= int(row[10])
		self.sbjct_start	= int(row[11])
		self.sbjct_end		= int(row[12])
		self.sbjct_len		= int(row[13])
		self.aln_query		= row[14]
		self.aln_sbjct		= row[15]
	
		
class ProjectFolders:
	"""folder structure of a project """
	
	def __init__(self, proj_home):

		self.home       = proj_home
		self.work       = "%s/work"        %  proj_home
		self.out        = "%s/output"      %  proj_home

		#working folders
		self.annotate   = "%s/annotate"    %  self.work
		self.lineage    = "%s/lineage"     %  self.work
		self.phylo      = "%s/phylo"       %  self.work
		self.internal   = "%s/internal"    %  self.work

		#second-level
		self.vgene      = "%s/vgene"       %  self.annotate
		self.jgene      = "%s/jgene"       %  self.annotate
		self.last       = "%s/last_round"  %  self.lineage
		self.beast      = "%s/beast"       %  self.phylo

		#output
		self.seq        = "%s/sequences"   %  self.out
		self.tables     = "%s/tables"      %  self.out
		self.plots      = "%s/plots"       %  self.out
		self.logs       = "%s/logs"        %  self.out
		self.rates      = "%s/rates"       %  self.out
		
		#second-level
		self.aa         = "%s/amino_acid"  %  self.seq
		self.nt         = "%s/nucleotide"  %  self.seq


#
# -- END -- class defination
#



#
# -- BEGIN -- argument process functions
#

def processParas(para_list, **keywords):
	""" 
	process the parameter information, all parameter start with "-paraname"
	return a dictionary (paraName, paraValue)
	remove the first parameter which is the program name
	"""
	
	para_list 		= para_list[1 :]
	kwgs, values 	= para_list[ :: 2], para_list[1 :: 2]

	if len(kwgs) != len(values):
		print "number of keywords and values does not equal"
		sys.exit(0)
	
	kwgs 	= map(lambda x : keywords[x[1 :]], kwgs)
	values 	= map(evalValues, values)

	#check for multiple uses of the same keyword and covert those values into a list
	for k in set(kwgs):
		if kwgs.count(k) > 1:
			idx     = [ i for i,x in enumerate(kwgs) if x==k ]
			valList = [ values[i] for i in idx ]
			values[idx.pop(0)] = valList #pop removes first occurence from idx list so we can only eliminate later occurences in next 2 lines
			kwgs    = [ x for i,x in enumerate(kwgs)   if i not in idx ]
			values  = [ x for i,x in enumerate(values) if i not in idx ]
			

	return dict(zip(kwgs,values))
	


def evalValues(v):
	"""Evaluate strings and return a value corresponding to its real type (int, float, list, tuple) """
	
	try:	return eval(v)
	except:	return v



def getParas(my_dict, *args):
	"""return the value of the argument """	
	
	if len(args) == 1:	return my_dict[args[0]]
	else:	return (my_dict.get(arg) for arg in args)
		
		
def getParasWithDefaults(my_dict, default_dict, *args):
	"""return the value of the argument or the default value if it wasn't specified on the command line"""	
	
	if len(args) == 1:	return my_dict.get(args[0],default_dict.get(args[0]))
	else:	return (my_dict.get(arg,default_dict.get(arg)) for arg in args)
		
		
#		
# -- END -- argument process functions
#


#
# -- BEGIN --  folder and file methods
#

def fullpath2last_folder(s):
	"""get immediate parent folder"""
	
	return s[s.rindex("/") + 1 :]


def create_folders(folder, force=False):

	old_wd = os.getcwd()
	os.chdir(folder)

	if (os.path.isdir("work")):
		if not force:
			sys.exit("Working directory already exists. Please use the -f(orce) option to re-intiate an analysis from scratch.\n")
		else:
			try:
				shutil.rmtree("work")
			except:
				sys.exit("Cannot remove old working directory, please delete manually and restart.\n")
			
	if (os.path.isdir("output")):
		if not force:
			sys.exit("Output directory already exists. Please use the -f(orce) option to re-intiate an analysis from scratch.\n")
		else:
			try:
				shutil.rmtree("output")
			except:
				sys.exit("Cannot remove old output directory, please delete manually and restart.\n")
			

	os.mkdir("work")
	os.mkdir("output")

		
	# Create working folders
	for subfolder in ALL_FOLDERS:
		try:
			os.mkdir(subfolder)
			
		except:		# may need to delete old folders
			print "FOLDER EXISTS: %s"%subfolder
			
	
	os.chdir(old_wd)
	return ProjectFolders(folder)




def check_folder(folder, sub_folder):
	"""return False if sub_folder does not exist; otherwise return full path of subfolder"""

	full_path = "%s/%s/" %(folder, sub_folder)
	
	if os.path.exists(full_path):
		return full_path

	return False


	

def get_files_format(folder, format = ""):
	"""return all files in specified folder with particular format; if format is not specified, return all files"""

	files = os.listdir(folder)
	
	if len(format) > 0:
		files = [x for x in files if x.endswith(format)]
	
	return sorted(files)




def get_files_format_fullpath(folder, format = ""):
	"""return all files in specified folder with particular format; if format is not specified, return all files"""

	files = os.listdir(folder)
	
	if len(format) > 0:
		files = [x for x in files if x.endswith(format)]
	
	files = ["%s/%s" %(folder, x) for x in files]
	
	return sorted(files)




def create_subfolder(folder, subfolder):
	"""create subfolder under current folder, nothing will be done if subfolder already exists"""
	
	full_path = "%s/%s/" %(folder, subfolder)
	print "CREATING FOLDER: %s" %full_path

	try:
		os.mkdir(full_path)

	except:		# normally the path exists, no action should be taken
		pass;
	
	return full_path


def full_path2file_name(full_path):
	ind = full_path.rindex("/")
	return full_path[ind + 1 :]


def parse_name(s):
	"""retrieve file name from full path """
	
	return s[s.rindex("/") + 1 : s.rindex(".")]




def check_files(paired_files):
	""" check whether fasta/qual/alignment file number and names agree with each other"""
	
	for values in paired_files:
		
		name_set = set()
		
		for v in values:
			name_set.add( parse_name(v) )
	
		if len(name_set) > 1 :
			print " ****** FATAL ERROR: FILE NAME DON'T MATCH ****** "
			print [parse_name(x) for x in values]
			print " ******  PLEASE CHECK YOUR FILES ******"
			sys.exit(0)


def create_array_job(target, name, resources, command, number, outfile=""):
	handle = open(target, "w")
	handle.write("#!/bin/bash\n")
	handle.write("#$ -t 1-%d\n" % number)
	handle.write("#$ -N %s\n" % name)
	handle.write("#$ -l %s\n" % resources)
	handle.write("#$ -cwd\n")
	if (output is not ""):
		handle.write("#$ -o %s\n" % output)
	handle.write("%s\n" % command)

	handle.close()
	os.system("qsub %s" % target)


def write_pbs_file(folder, library):
	
	#adjust blast options for V vs D/J
	for infile in glob.iglob("%s/*.fasta"%folder):
		
		fhead, ftail = os.path.splitext(infile)
		f_ind = fhead[ fhead.rindex("_") + 1 :]
		
		# GERMLINE pbs
		germ_file 		= open("%s/blast%s.sh" %(folder, f_ind), "w")
		job_name 		= "g%s" %f_ind
		outfile 		= "%s/%s.txt" %(folder, fhead)
		infile_fullpath 	= "%s/%s/%s" %(os.getcwd(), folder, infile)
		pbs_string 		= PBS_STRING %(job_name, outfile, options, library, infile_fullpath)
		
		germ_file.write(pbs_string)
		germ_file.close()
		
		
def write_auto_submit_file(folder_tree):
	"""write shell script to automatically submit all jobs to cluster"""
	
	print "generating auto submit file..."
	
	handle 		= open("%s/submit_all.sh" %folder_tree.jobs, "w")
	infiles 	= get_files_format(folder_tree.pbs, "sh")
	
	for infile in infiles:
		handle.write("qsub %s/%s\n" %(folder_tree.pbs, infile))
		#handle.write("sleep 2\n")
		handle.write("\n")



def write_read_length(folder_tree, l, total): # l is an array
	"""record read distribution"""
	
	prj_name 	= folder_tree.home[folder_tree.home.rindex("/") :]
	outfile 	= "%s/%slength_distribution.txt" %(folder_tree.data, prj_name)
	writer 		= csv.writer(open(outfile, "w"), delimiter = sep)
	
	writer.writerow(["length", "count", "percentile"])
	
	max_len = len(l)
	for ind, rlen in enumerate(l):
		aline = [ind, rlen, float(rlen) / total * 100]

	
def load_divid(f):
	reader, result = csv.reader(open(f, "rU"), delimiter = sep), dict()
	reader.next()
	for row in reader:
		result[row[0]] = row[ 1 :]
	return result
	
	
def load_primer(f):
	reader, result = csv.reader(open(f, "rU"), delimiter = sep), dict()
	reader.next()
	for row in reader:
		result[row[0]] = row[1]
	return result	
		

#
# -- END -- folder and file methods 
#


#
# -- BEGIN -- qual file methods
#

def generate_phred_scores(f):
	"""parse quality file and generate phred scores one read a time"""
	
	handle = open(f, "rU")
	old_id, old_quals = get_454_fasta_id(handle.next()), []

	for line in handle:
		if line.startswith(">"):					# start of a new entry
			new_id = get_454_fasta_id(line)

			yield MyQual(old_id, old_quals)
			old_id, old_quals = new_id, []
			
		else:		# quality score lines	
			line 		= line[ : -1].strip().split(" ")	# remove "\n", remove " " at both ends, split by " "
			line 		= [x for x in line if x != ""]
			old_quals 	+= map(int, line)
			
	yield MyQual(old_id, old_quals)
	

def generate_quals_folder(folder):
	qual_files = glob.glob("%s/*.qual" %folder)
	for qual_file in qual_files:
		for myqual in generate_phred_scores(qual_file):
			yield myqual
	


def load_quals(f):
	"""reutrn sequence ID and quality scores in dicitonary format"""

	result = dict()
	
	for seq_id, quals, seq_len in generate_phred_scores(f):
		result[seq_id] = quals
		
	return result
		

#
# -- END -- qual file methods
#


#
# -- BEGIN -- FASTA file and sequence methods
#

def has_pat(s, pat):
	has, start, end = False, -1, -1
	matches = re.finditer(pat, s)
	for match in matches:
		has, start, end = True, match.start(), match.end()
	return has, start, end
	

def write_seq2file(myseq1, myseq2, f):
	""" write two MySseq objects to given file"""
	try:
		handle = open(f, "w")
		handle.write(">%s\n" %myseq1.seq_id)
		handle.write("%s\n\n" %myseq1.seq)
	
		handle.write(">%s\n" %myseq2.seq_id)
		handle.write("%s\n\n" %myseq2.seq)
		handle.close()
	except:
		print type(myseq1)
		print type(myseq2)



def load_seqs_in_dict(f, ids):
	"""
	load all sequences in file f is their id is in ids (list or set or dictionary)
	"""
	result = dict()
	for entry in SeqIO.parse(open(f, "rU"), "fasta"):
		if entry.id in ids:
			result[entry.id] = entry
			
	return result
	


def generate_read_in_dict(f, d):
	print "loading sequence from %s..." %f
	reader, total = SeqIO.parse(open(f, "rU"), "fasta"), 0
	for entry in reader:
		if entry.id in d:
			total += 1
			if total % 10 ** 4 == 0:
				print "%d..." %total

			yield MySeq(entry.id, entry.seq)
			
			
	print "total %d loaded..." %total

#def get_202_germline():
#	for germ in SeqIO.parse(open(GERM_HEAVY, "rU"), "fasta"):
#		if germ.id == "IGHV1-2*02":
#			return MySeq("V202", germ.seq)
		
def get_germline(source, is_light):
	for germ in SeqIO.parse(open(dict_germ_db[is_light], "rU"), "fasta"):
		if germ.id == source:
			return MySeq(source, germ.seq)
		
def load_germlines():
	"""return germline gene ID and sequence in a dictionary """

	print "loading germline sequences..."
	
	reader, result = SeqIO.parse(open(GERM_DB, "rU"), "fasta"), dict()

	for entry in reader:
		result[entry.id] = MySeq(entry.id, entry.seq)
	
	return result
	

def load_fasta_dict(f):
	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), dict()
	for entry in reader:
		result[entry.id] = entry
	return result


def load_fastas_in_list(f, l):
	
	print "loading reads from %s as in given list..." %f
	reader, result, good = SeqIO.parse(open(f, "rU"), "fasta"), dict(), 0

	for entry in reader:
		if entry.id in l:

			#changed to match load from set CAS 20121004
			result[entry.id] = entry
			good += 1
			if good == len(l): break

			#myseq = MySeq(entry.id, entry.seq)
			#myseq.desc = entry.description
			#result[entry.id] = myseq

	print "%d loaded...." %len(result)
	return result
	
def load_fastas_in_set(f, s):
	print "loading reads from %s as in given set..." %f
	reader, dict_reads, good = SeqIO.parse(open(f, "rU"), "fasta"), dict(), 0
	for ind, entry in enumerate(reader):
		if (entry.id in s or int(entry.id) in s):
			dict_reads[entry.id] = entry
			good += 1
			
			if good == len(s):
				break
	print "%d loaded..." %len(dict_reads)
	return dict_reads


def load_list_fastas(f, l):
	print "loading fasta from %s... " %f
	#print " ".join(l)
	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), []
	for entry in reader:
		if entry.id in l:
			myseq = MySeq(entry.id, entry.seq)
			myseq.desc = entry.description
			result.append(myseq)
	return result


def load_fastas_with_Vgene(f, v):
	print "loading reads from %s assigned to %s..." %(f,v)
	reader, dict_reads = SeqIO.parse(open(f, "rU"), "fasta"), dict()
	for entry in reader:
		if re.search(v, entry.description):
			dict_reads[entry.id] = entry

	print "%d loaded..." %len(dict_reads)
	return dict_reads


def load_heavy_fasta_V(f, l):
	print "loading fasta from %s... " %f
	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), []
	for entry in reader:
		if entry.id in l:
			# truncate V gene
			has_wgxg, wgxg_start, wgxg_end = has_pat(str(entry.seq), pat_nuc_wgxg)
			has_cxrk, cxrk_start, cxrk_end = has_pat(str(entry.seq), pat_nuc_cxrk)
			entry.seq = entry.seq[ : cxrk_end]
			
			myseq = MySeq(entry.id, entry.seq)
			myseq.desc = entry.description
			result.append(myseq)
	return result

	
def load_fastas(f):
	"""return gene ID and sequences in a dictionary"""
	print "loading sequence info from %s..." %f

	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), dict()

	for entry in reader:
		#myseq = MySeq(entry.id, entry.seq)
		#myseq.desc = entry.description
		result[entry.id] = entry
	
	return result

def load_fastas_translation(f):
	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), []

	for entry in reader:
		myseq = MySeq(entry.id, entry.seq.translate())
		result.append(myseq)
	return result

def fix_N(seq):
	# remove "N" at the beginning
	while seq.startswith("N"):
		seq = seq[ 1 : ]

	# remomve "N" at the end
	while seq.endswith("N"):
		seq = seq[ : -1]
	
	# remove "N" in the middle
	while seq.find("N") > 0:
		pos = seq.index("N")
		seq = "".join([seq[ : pos], seq[pos - 1], seq[pos + 1 :]])

	return seq
	

def generate_read_fasta(f):
	"""read fasta file and yield one reads per time """
	
	reader = SeqIO.parse(open(f, "rU"), "fasta")
	for entry in reader:
		yield entry


def generate_read_fasta_folder(fastas, type=1):

	for fasta_file in fastas:

		filetype = "fasta"
		if re.search("\.(fq|fastq)$", fasta_file) is not None:
			filetype = "fastq"

		for entry in SeqIO.parse(open(fasta_file, "rU"), filetype):

			yield MySeq(entry.id, entry.seq), None, fasta_file

			#if we reimplement qual handling uncomment next section
			#if filetype == "fastq":
			#	yield MySeq(entry.id, entry.seq, desc), MyQual(entry.id, entry.letter_annotations["phred_quality"]), fasta_file	
			

def generate_reads(f):
	"""read fasta file and yield one reads per time """

	reader = SeqIO.parse(open(f, "rU"), "fasta")
	for entry in reader:
		# entry.seq = fix_N(entry.seq.tostring().upper())
		myseq = MySeq(entry.id, entry.seq)
		
		try:
			desc = entry.description
			myseq.desc = desc
		except:
			pass;
		
		yield myseq

def get_fasta_len_str(f):
	print "loading fasta length from file: %s..." %(f), 
	result = dict()
	for entry in SeqIO.parse(open(f, "rU"), "fasta"):
		result[entry.id] = len(entry.seq)
	print "done! %d loaded..." %len(result)
	return result

def str_reverse_complement(s):
	"""reverse complement a sequence in string format """
	
	return str(Seq.Seq(s).reverse_complement())
	
	

def get_fasta_len(f):
	""" return sequence id and length in dictionary format  """
	print "reading %s... for read length" %f
	
	reader, result = SeqIO.parse(open(f, "rU"), "fasta"), dict()
	
	for entry in reader:
		result[int(entry.id)] = len(entry.seq)
		
	print "total %d reads/sequences in %s." %(len(result), f)
	
	return result
	
	
	
def check_fasta_qual_pair(fa_name, qu_name):
	return fa_name[ : fa_name.rindex(".")] == qu_name[ : qu_name.rindex(".")]


"""		non-compatible izip on Titan server
def generate_fasta_qual_pair(folder):
	
	
	fastas, quals = get_files_format_fullpath(folder, "fasta"), get_files_format_fullpath(folder, "qual")	
	pairs = zip(fastas, quals)
	
	for fasta_file, qual_file in pairs:
		print "fasta_file: %s \r\nqual_file: %s" %(fasta_file, qual_file)
		if not check_fasta_qual_pair(fasta_file, qual_file):
			print "FATAL ERROR: NAME DOES NOT MATCH"
			print fasta_file, qual_file
			#sys.exit(0)
			
		
		for myseq, myqual in izip(generate_read_fasta(fasta_file), generate_phred_scores(qual_file)):
			if (myseq.seq_id == myqual.qual_id) and (myseq.seq_len == myqual.qual_len):
				#yield myseq, myqual
				pass;
				
			else:
				print "FATAL ERROR: READ: ID %s LEN: %d DOES NOT MATCH QUAL: ID %s LEN %d" %(myseq.seq_id, myseq.seq_len, myqual.qual_id, myqual.qual_len)
				#sys.exit(0)
			
			yield myseq, myqual
"""

def cut_from_stop_codon(s):
	if s.find("*") >= 0:
		return s[ : s.index("*")]
	return s



def translate_a_sequence(s):
	"""translate nucleotide to protein"""  # this method is redundant, but we wannt deal with "N"
	s = s.upper()
	if s.find("N") > 0:
		return None
	codons = ["".join(x) for x in zip(s[ :: 3], s[1 :: 3], s[2 :: 3])]

	return "".join([dict_codon2aa[x][0] for x in codons])
	


def get_a_seq(f, seq_id):
	reader = SeqIO.parse(open(f, "rU"), "fasta")
	for entry in reader:
		if entry.id == seq_id:
			return MySeq(entry.id, entry.seq)
	
	print "No sequence: %s" %seq_id
	sys.exit(0)



def get_element(sequence, pattern, min_start):
	
	for match in re.finditer(pattern, sequence):
		start, end = match.start(), match.end()			# , match.group(0)
		if start > min_start:
			return start, end

	return -1, -1		# motif not found
	


def find_all_elements(sequence, pattern):
	result = []
	for match in re.finditer(pattern, sequence):
		start, end = match.start(), match.end()			# , match.group(0)
		result.append((start, end))
	return result
	


def translate_partial_seq(seq, end):
	region_len 		= end / 3 * 3
	start 			= end - region_len
	coding_region 	= seq[start : end].upper()
	
	# ambiguous nucleotide
	if coding_region.find("N") > 0:		
		return None	
	

	# a good translation
	protein = translate_a_sequence(coding_region)
	if protein.find("*") < 0 :
		return protein

	# bad translation but not curatable
	if not coding_region.endswith("TTT"):
		return None


	# try to remove a nucleotide OVERCALL
	end 			-= 1
	region_len 		= end / 3 * 3
	start 			= end - region_len
	coding_region 	= seq[start : end]
	protein 		= translate_a_sequence(coding_region)

	if protein.find("*") < 0:
		return protein

	# not OVERCALL, try UNDERCALL
	end 			+= 2
	region_len 		= end / 3 * 3
	start 			= end - region_len
	coding_region 	= seq[start : end]
	protein 		= translate_a_sequence(coding_region)
	
	if protein.find("*") < 0:
		# Sorry, this cannot be fixed by OVER-/UNDER-CALL
		return protein
	
	return None				
	


def curate_wgxg(s):
	middle, result 	= s[3 : -3], None		# W..G
	mid_len			= len(middle)
	
	if not middle.startswith("GG"):		# "G"  undercall at the begging
		middle = "G" + middle
	
	if len(middle) > 6 and middle.endswith("G"):
		middle = middle[ : -1]
		
	if len(middle) > 6:
		return None
	
	return "".join(["W", "G", dict_codon2aa[middle[ -3 :]][0], "G"])
	



def load_single_fasta(f, gid):
	""" retrieve single gene info from fasta file """
	for entry in SeqIO.parse(open(f, "rU"), "fasta"):
		if entry.id == gid:
			return MySeq(entry.id, entry.seq)
			
	print "Gene: %s not found in file: %s" %(gid, f)

def load_single_germ(ab_id):
	for ab in SeqIO.parse(open(GERM_DB, "rU"), "fasta"):
		if ab.id == ab_id:
			return ab
	return None
	
def load_single_germ_db(GERM_DB, ab_id):
	for ab in SeqIO.parse(open(GERM_DB, "rU"), "fasta"):
		if ab.id == ab_id:
			return ab
	return None

def load_single_germ_light(abid):
	for Ab in SeqIO.parse(open(GERM_LIGHT, "rU"), "fasta"):
		#print Ab.id
		if Ab.id == abid:
			return Ab
	return None

def load_single_native(ab_id):
	for ab in SeqIO.parse(open(NAT_DB, "rU"), "fasta"):
		if ab.id == ab_id:
			return ab
	return None
	
def load_single_native_light(abid):
	for Ab in SeqIO.parse(open(NAT_LIGHT, "rU"), "fasta"):
		if Ab.id == abid:
			return Ab
	return None

def format_seq(s):
	s_len = len(s)
	starts = range(0, s_len, 60)
	return "\n".join([s[start : start + 60] for start in starts])	

def load_native_cdrh3_dict():
	print "loading native CDRH3s ..."
	result = dict()
	for entry in SeqIO.parse(open(NAT_H_CDR3, "rU"), "fasta"):
		#if entry.id in l:
		result[entry.id] = entry
	print "%d native CDRH3s loaded" %len(result)
	return result
	
	
def load_list_native_cdrh3_dict(l):
	print "loading native CDRH3s ... %s" %",".join(l)
	result = dict()
	for entry in SeqIO.parse(open(NAT_H_CDR3, "rU"), "fasta"):
		if entry.id in l:
			result[entry.id] = entry
	print "%d native CDRH3s loaded" %len(result)
	return result

def load_native_in_dict():
	result = dict()
	for ab in SeqIO.parse(open(NAT_HEAVY, "rU"), "fasta"):
		result[ab.id] = ab
	return result

def load_native_in_dict2(is_light):
	result = dict()
	for ab in SeqIO.parse(open(dict_nat_db[is_light], "rU"), "fasta"):
		result[ab.id] = ab
	return result


def get_subfolders(f):
	return sorted([x for x in os.listdir(f) if os.path.isdir(x)])
	
#
# -- END -- FASTA file and sequence methods 
#



#
# -- BEGIN -- alignment file methods 
#

# clustalw

def remove_start_insertion(tst_seq, ref_seq):
	while ref_seq.startswith("-"):
		tst_seq = tst_seq[ 1 :]
		ref_seq = ref_seq[ 1 :]
	return tst_seq, ref_seq


def remove_an_insertion(tst_seq, ref_seq):
	ind = ref_seq.index("-")
	tst_seq = tst_seq[ : ind] + tst_seq[ind + 1 :]
	ref_seq = ref_seq[ : ind] + ref_seq[ind + 1 :]
	return tst_seq, ref_seq


def remove_dashes(s):
	s = s.split("-")
	return "".join(s)

def remove_pair_insertion_ref(tst_seq, ref_seq):
	while ref_seq.find("-") >= 0:
		ind = ref_seq.index("-")
		tst_seq, ref_seq = tst_seq[ : ind] + tst_seq[ ind + 1 :], ref_seq[ : ind] + ref_seq[ ind + 1 :]
	return tst_seq, ref_seq

def remove_end_insertion(tst_seq, ref_seq):
	while ref_seq.endswith("-"):
		tst_seq, ref_seq = tst_seq[ : -1], ref_seq[ : -1]
	return tst_seq, ref_seq


def removeMatchedGaps(ref, test):
	#by CAS 2012-12-05
	gaps = [m.start() for m in re.finditer("-", ref)]
	for dash in reversed(gaps):
		if dash < len(test):
			if (test[dash] == "-"):
				test, ref = test[ : dash] + test[ dash+1 :], ref[ : dash] + ref[ dash+1 :]	
	return ref,test


def trimAllEndGaps(first, second):
	#by CAS 2012-12-05
	#added iterative loops 2012-12-12
	#     (may not be worth the overhead for rare cases?)

	while any( x is not None for x in [re.match("-+", first), re.match("-+", second)] ):
		m = re.match("-+", first)
		if m:
			first,second = first[m.end():], second[m.end():]
		n = re.match("-+", second)
		if n:
			first,second = first[n.end():], second[n.end():]

	while any( x is not None for x in [re.search("-+$", first), re.search("-+$", second)] ):
		m = re.search("-+$", first)
		if m:
			first,second = first[:m.start()], second[:m.start()]
		n = re.search("-+$", second)
		if n:
			first,second = first[:n.start()], second[:n.start()]

	return first, second		


def prepare_clustal_jobs(prj_tree, germ):
	"""
	prepare pbs files for each of the fasta files in clustal_fasta
	"""
	
	# load all fasta files in clustal_fasta folder
	infiles = glob.glob("%s/*.fa" %prj_tree.clustal_fasta)
	for infile in infiles:
		head, tail 	= os.path.splitext(infile)
		f_ind 		= head.split("_")[ -1 ]
		
		for native in DICT_H_NATIVES[germ]:
			handle = open("%s/%s_%s.sh" %(prj_tree.clustal_pbs, native, f_ind), "w")
			handle.write("#!/bin/bash\n")
			handle.write("#$ -N %s_%s\n" %(native, f_ind))
			handle.write("#$ -l mem=100M,time=10:00:00\n")
			handle.write("#$ -cwd\n")
			handle.write("#$ -o %s/%s_%s.txt\n\n" %(prj_tree.clustal_data, native, f_ind))
			
			handle.write("1.24-do_clustal.py -g %s -n %s -i %s\n" %(germ, native, infile))
			handle.close()

def is_valid_native_aln(align_file):
	try:
		alignment 	= AlignIO.read(align_file, "clustal")
	except:				# implement for Titan server with BioPython 1.52
		alignment 	= AlignIO.read(open(align_file, "rU"), "clustal")
	read 	= str(alignment[0].seq)		# first sequence is test
	total_deletion, ind = 0, 0
	while True:
		if read[ind] == "-":
			total_deletion += 1
			ind += 1
		else:
			break
	ind = 1
	while True:
		if read[-ind] == "-":
			total_deletion += 1
			ind += 1
		else:
			break
	#print "total deletion: %d" %total_deletion
	return total_deletion < MAX_NATIVE_CLUSTAL_DELETION
	
def get_clustal_matches(align_file):
	try:
		alignment 	= AlignIO.read(align_file, "clustal")
	except:				# implement for Titan server with BioPython 1.52
		alignment 	= AlignIO.read(open(align_file, "rU"), "clustal")
	seq1, seq2 	= str(alignment[0].seq), str(alignment[1].seq)
	zip_seqs 	= zip(seq1, seq2)
	zip_seqs	= [(x, y) for x, y in zip_seqs if x !="-" and y != "-"]		# remove insertions from either one
	matches		= sum([x == y for x, y in zip_seqs])						# count total matches
	#identity 	= round(float(matches) / my_ref_len * 100, 2)
	
	mismatches	= sum([x != y for x, y in zip_seqs])
	#divergence	= round(float(mismatches) / my_ref_len * 100, 2)
	
	return matches, mismatches
		
	

def parse_pair_clustal(align_file, my_ref_len, get_coverage=0):
	""" parse paired alignment file """
	try:
		alignment 	= AlignIO.read(align_file, "clustal")
	except:				# implement for Titan server with BioPython 1.52
		alignment 	= AlignIO.read(open(align_file, "rU"), "clustal")
	seq1, seq2 	= str(alignment[0].seq), str(alignment[1].seq)
	zip_seqs 	= zip(seq1, seq2)
	zip_seqs	= [(x, y) for x, y in zip_seqs if x !="-" and y != "-"]		# remove insertions from either one
	matches		= sum([x == y for x, y in zip_seqs])						# count total matches
	identity 	= round(float(matches) / my_ref_len * 100, 2)
	
	mismatches	= sum([x != y for x, y in zip_seqs])
	divergence	= round(float(mismatches) / my_ref_len * 100, 2)

	if get_coverage > 0:
		i =len(filter(lambda let: let == "-", seq1))
		j =len(filter(lambda let: let == "-", seq2))
		my_ref_len = min(len(seq1)-i,len(seq2)-j) #only require one-sided coverage (catches in-dels), so use shorter length
		identity = float(matches) / my_ref_len * 100
		coverage = float(len(zip_seqs)) / my_ref_len * 100
		return identity, coverage
	else:
		return identity, divergence


def do_clustalw(ref, tst, fa_file=""):

	fhead=fa_file
	if fa_file == "":
		fhead = generate_random_string(8)
		fa_file = "%s.fa" % fhead

	write_seq2file(ref, tst, fa_file)
	
	#update to choose method
	cline = ClustalwCommandline(clustalw, infile = fa_file)
	cline()

	return fhead


def get_clustalw2_score(status):
	status = status[-2].split("\n")
	
	score1 = status[12].split(" ")[-1]
	score2 = status[20].split(" ")[-1]
	return score1, score2
	
	for item in status:
		if item.find("Alignment Score") >= 0:
			return item.split(" ")[-1]	
						
		

def do_clustalw_v2(tst, ref, fa_file):
	#print "write to file...."
	write_seq2file(tst, ref, fa_file)

	#print "do clustal..."
	# Titan cluster
	if HOSTNAME.find("titan") >= 0:
		#print "titan", fa_file
		#print "clustalw2 %s" %fa_file
		commands.getstatusoutput("clustalw2 %s" %fa_file)

	else:
		cline = ClustalwCommandline(clustalw, infile = fa_file)
		status = cline()
		return get_clustalw2_score(status)
	#print "clustal is done!"
	
def do_clustalw2_v3(ref, tst, fa_file):
	#print "inside function..."
	handle = open(fa_file, "w")
	SeqIO.write(ref, handle, "fasta")
	SeqIO.write(tst, handle, "fasta")
	handle.close()
	#print "file written"
	if HOSTNAME.find("titan") >= 0:
		#print "titan", fa_file
		#print " %s" %fa_file
		commands.getstatusoutput("clustalw2 %s" %fa_file)

	else:
		try:
			cline = ClustalwCommandline(clustalw, infile = fa_file)
			status = cline()
			#print status
		except:
			print "clustalw2 wrong!"
			sys.exit(0)
		
	return

def write_file_2seq_align(ref, tst, f):
	handle = open(f, "w")
	SeqIO.write(ref, handle, "fasta")
	SeqIO.write(tst, handle, "fasta")
	handle.close()

def do_2seq_clustalw2_alignment(ref, tst, fa_file):
	write_file_2seq_align(ref, tst, fa_file)
	if HOSTNAME.find("titan") >= 0:
		#print "titan", fa_file
		#print "clustalw2 %s" %fa_file
		commands.getstatusoutput("clustalw2 %s" %fa_file)

	else:
		cline = ClustalwCommandline(clustalw, infile = fa_file)
		status = cline()

def get_clustal_align_seq_score(aln_file):
	alignment 	= AlignIO.read(aln_file, "clustal")
	seq1, seq2 	= str(alignment[0].seq), str(alignment[1].seq)
	seqz = zip(seq1, seq2)
	total_id = sum([ x == y for (x, y) in seqz])
	percent = float(total_id) / len(seq1) * 100
	return seq1, seq2, percent
	
	
def score_alignment(aln1, aln2):
	seqz = zip(aln1, aln2)
        total_id = sum([ x == y and x != "-" for (x, y) in seqz])
	align_len = min(len("".join(aln1.split("-"))),len("".join(aln2.split("-"))))
	percent = float(total_id) / align_len * 100
        return percent


def chaims_alignment_score(one, two):
	seqz=zip(one.upper(),two.upper())
	matched = sum([ x == y for (x, y) in seqz])
	ungapped = re.split("-+",one) + re.split("-+",two)
	align_len = len("".join(ungapped))
	num_gaps = len(ungapped) - 2
	gap_len = len(re.findall('-',one) + re.findall('-',two))

	#print [matched, num_gaps, gap_len, align_len, len("".join(re.split("-+",one))), len("".join(re.split("-+",two)))]

	#score = float((matched * 2) - (5 * num_gaps))
	#score = float((matched * 2) - (10 * num_gaps)) / (align_len ** 0.5)
	score = float((matched * 2) - (10 * num_gaps)) / align_len
	
	return score

	
def remove_ref_insertion_both_ends(r, t):
	while r.startswith("-"):
		r, t = r[ 1 :], t[ 1 :]
	while r.endswith("-"):
		r, t = r[ : -1], t[ : -1]	

	return r, t


def remove_ref_insertion(ref, tst):
	while ref.find("-") > 0:
		index = ref.index("-")
		ref = ref[ :index] + ref[index + 1 :]
		tst = tst[ :index] + tst[index + 1 :]
	return ref, tst


def parse_dnd(f):
	node_list, node_stack, dict_node_children, node_ind = [], [], dict(), 0
	def add_child(node):
		for parent in node_stack:
			dict_node_children[parent].add(node)

	for ind, line in enumerate(open(f, "rU")):
		line = line.strip()
		print "line:", ind, line
		print "stack:", node_stack

		if line.startswith("("):
			# this is an new intermediate node
			new_node = "N%07d" %node_ind
			# record current intermediate node
			node_list.insert(0, new_node)
			# add current node as child of previous nodes
			#add_child(new_node)
			# push current node into stack
			node_stack.insert(0, new_node)
			# initiate current node's children
			dict_node_children[new_node] = set()
			node_ind += 1


		elif line.startswith(":"):
			# the branch length of previous node
			pass;
		elif line.find(":") > 0:
			# this is a leaf node
			leaf = line.split(":")[0]
			# add this leaf to all current intermediate nodes
			add_child(leaf)
		if line.endswith(")") or line.endswith(";"):		# ";" indicates the last line
			# end of current intermediate node
			node_stack.pop(0)

		#print "node and children:"
		#print "-------------"	
		#for node in node_list:
		#	print node, dict_node_children[node]
		#print "-------------"
		#print

	# print "Node Stack is empty:", len(node_stack) == 0, node_stack
	leafs_list = dict_node_children.values()
	len_leafs  = sorted([(len(x), x) for x in leafs_list])
	sorted_leafs = [leaf for (len_leaf, leaf) in len_leafs]
	return sorted_leafs


def get_germ_vague(f, germ):
	print "loading germs of %s from file: %s..." %(germ, f)
	reader, result = csv.reader(open(f, "rU"), delimiter = sep), set()
	reader.next()
	for row in reader:
		read_id, germ_id = int(row[0]), row[1].strip()
		if germ in germ_id:
			result.add("%08d" %read_id)
	return result
	

def get_germ_exact(f, germ):
	print "loading germs of %s from file: %s..." %(germ, f)
	reader, result = csv.reader(open(f, "rU"), delimiter = sep), set()
	header = reader.next()
	for row in reader:
		read_id, germ_id = int(row[0]), row[1].strip()
		if germ == germ_id:
			result.add("%08d" %read_id)
	return result
	


def generate_top_hits(f):
	"""return top hit of the alignment """		#****** This method will be obsolete soon ******#
	
	reader, old_id = csv.reader(open(f, "rU"), delimiter = sep), ""
	
	for row in reader:
		qid = row[0]
		
		if qid != old_id:
			sid, qstart, qend, sstart, send, strand, aln_len, mismatch, gap = row[1], int(row[6]), int(row[7]), int(row[8]), int(row[9]), "+", int(row[3]), int(row[4]), int(row[5])
		
			if (qstart - qend) * (sstart - send) < 0:
				strand 			= "-"
				sstart, send 	= send, sstart
			
			yield qid, sid, strand, qstart, qend, sstart, send, aln_len, mismatch, gap, row
			old_id = qid
	
		else:
			pass;


# BLAST

def get_col_id(headers, header):
	for ind, h in enumerate(headers):
		if header == h:
			return ind
	return -1


def load_divid_native_sorted(f, nat, mg, ng, mn, nn):
	"""
	"nat" is a column header in file "f", where the read to native identities are stored.
	mg, ng: mx_germ_divergence, mn_germ_divergence
	mn, nn: mx_nat_identity, min_nat_identity
	Load all reads with native identity in [nn : mn] and germ divergence in [ng : mg]
	Sort them based on native identity from large to small and germ divergence from large to small
	return sorted read IDs 
	"""
	reader = csv.reader(open(f, "rU"), delimiter = sep)
	headers = reader.next()
	
	col_id = get_col_id(headers, nat)
	if col_id == -1:
		print "native antibody %s not in divergence/identity file" %nat
		sys.exit(0)
		
	divid_list, dict_read_divid = [], dict()
	for row in reader:
		read, germ_div, nat_id = row[0], float(row[1]), float(row[col_id])
		if ng <= germ_div <= mg and nn <= nat_id <= mn:
			divid_list.append((nat_id, germ_div, read))
			dict_read_divid[read] = (germ_div, nat_id)

	print "%d loaded... sorting" %len(divid_list)
	divid_list 		= sorted(divid_list)[ :: -1]
	sorted_reads 	= [x for (y, z, x) in divid_list]

	return sorted_reads, dict_read_divid
	


def get_germ_in_dict(f):
	# load germline assignment into dict
	print "loading germlines from file: %s..." %f
	result = dict()
	reader = csv.reader(open(f, "rU"), delimiter = sep)
	reader.next()
	for ind, row in enumerate(reader):
		read_id, germ, strand, aln_len = int(row[0]), row[1].strip(), row[12].strip(), int(row[3])
		if aln_len > 150:
			result[read_id] = (germ, strand)
		if ind % 10 ** 4 == 0:
			print "%d loaded..." %ind
	print "total: %d." %(ind + 1)
	return result

			
def gernerate_germ_assign_folder_v2(folder, germ, min_cov=.98, max_start = 3):
	"""
	return reads that was assigned to germline according to minimum coverage and maximum start
	"""
	infiles, total, min_len = glob.glob("%s/*.txt" %folder), 0, germ.seq_len * min_cov
	for infile in infiles:
		print "processing %s...." %infile
		reader = csv.reader(open(infile, "rU"), delimiter = sep)
		old_id, current_total = "", 0
		for row in reader:
			try:
				qid, sid, sstart, send, evalue = row[0], row[1], int(row[8]), int(row[9]), float(row[10])
				scov = abs(sstart - send) + 1
				
				if sid == germ.seq_id and sstart <= max_start and scov >= min_len:
					strand = "+"
					if sstart > send:
						strand = "-"
					
					yield qid, evalue, strand
					total += 1			
				
			except:
				pass
		print "total: %d" %total


def generate_germ_assign_file(f, germ):
	
	reader = csv.reader(open(f, "rU"), delimiter = sep)
	old_id, total = "", 0
	for row in reader:
		try:
			qid, sid, sstart, send, evalue = row[0], row[1], int(row[8]), int(row[9]), float(row[10])
			scov = abs(sstart - send) + 1
			
			if qid == old_id:
				pass
			elif sid == germ and evalue <= GERM_EVALUE and scov >= GERM_MINCOV:
				strand = "+"
				if sstart > send:
					strand = "-"
				
				old_id = qid
				yield qid, strand
				total += 1			
			
		except:
			pass
	print "total: %d" %total
			

def generate_germ_assign_folder(foler, germ, e_value, top=3):
	"""
	return reads assigned to germ with maximum evalue given
	"""
	infiles, total = glob.glob("%s/*.txt" %foler), 0
	for infile in infiles:
		print "processing %s...." %infile
		reader = csv.reader(open(infile, "rU"), delimiter = sep)
		old_id, current_total = "", 0
		for row in reader:
			#for ind, item in enumerate(row):
			#	print ind, item
				
			#sys.exit(0)
			try:
				qid, sid, sstart, send, evalue = row[0], row[1], int(row[8]), int(row[9]), float(row[10])
				if qid != old_id:
					old_id, current_total = qid, 1
				else:
					current_total += 1
			
				if current_total <= top and sid == germ and evalue <= e_value:
					strand = "+"
					if sstart > send:
						strand = "-"
				
				
					yield qid, evalue, strand
					total += 1
			except:
				pass;
		print "total: %d" %total

			

def get_dict_blast_tophit_verbose(f):
	"""return dictionay of verbose blast result"""

	print "retrieving top hits from file: %s..." %f
	reader, result, old_id = csv.reader(open(f, "rU"), delimiter = sep), dict(), ""
	
	reader.next()
	for row in reader:
		query_id = row[0]
		
		if query_id != old_id:
			my_alignment = MyAlignmentVerbose(row)
			result[query_id] = my_alignment
			
			old_id = query_id
	
	if query_id not in result:
		result[query_id] = my_alignment
	
	return result


def cut_five_end_blast(str_seq, my_aln):
	start, end = my_aln.qstart, my_aln.qend
	
	if my_aln.strand == "+":
		return str_seq[start - 1 :]
	elif my_aln.strand == "-":
		return str_seq[ : end].reverse_complement()


def template_replace_N(tst, ref):
	"""
	Replace n in tst with nucleotide on ref
	"""	
	zseq, result_tst, result_ref = zip(tst, ref), [], []
	for x, y in zseq:
		if x == "n" and y == "-":
			pass;				# skip ambiguous nucleotide
		else:
			if x == "n":
				x = y			# replace reads nucleotide with reference nucleotide
			
			# record current nucleotides in ref and tst
			result_tst.append(x)
			result_ref.append(y)
	return "".join(result_tst), "".join(result_ref)
	
	



def parse_blast_strand(blast_strand):
	# ('Plus', 'Plus') = "+"; ('Plus', 'Minus') = "-"
	if blast_strand[0] == blast_strand[1]:
		return "+"
	else:
		return "-"			
			
			

def parse_blast_gaps(blast_gaps):
	# (6, 304) = "6/304"
	return "%d/%d" %(blast_gaps[0], blast_gaps[1])


def replace_not_aligned(read, ref):
	z, result = zip(read, ref), []
	
	for x, y in z:
		if x in "N-":
			x = y
		result.append(x)

	return "".join(result) 

def replace_N(ref, germ, tst):
	if len(ref) != len(tst):
		print "not equal length: put: %d; read: %d" %(len(ref), len(tst))
		
	z, result = zip(ref, germ, tst), []	
	for r, g, t in z:
		if t == "N":
			if r == "N":
				t = g
			else:
				t = r
		
		result.append(t)
	
	return "".join(result)
	

def remove_dash(s):
	return "".join([x for x in s.split("-") if x != "-"])


def generate_hits_identity_diversity(f):
	""" 
	deal with processed top hits, recalcualte identity
	"""
	print "retrieving top hit information from %s..." %f
	reader, total = csv.reader(open(f, "rU"), delimiter = sep), 0
	reader.next()
	
	for row in reader:
		my_alignment = MyAlignment(row)
		
		my_alignment.qlen 		= int(row[-2])
		my_alignment.slen 		= int(row[-1])
		my_alignment.strand 	= row[-3]
		
		matches = my_alignment.alignment - my_alignment.mismatches - my_alignment.gaps
		
		real_identity 	= round(float(matches) / my_alignment.slen * 100, 2)
		divergence		= round(100 - real_identity, 2)
		
		my_alignment.real_id 		= real_identity
		my_alignment.divergence 	= divergence
		
		total += 1
		
		if total % 10 ** 4 == 0:
			print "%d processed..." %total
		
		yield my_alignment


"""	
def generate_array(dict_div_id):
	my_array, total = zeros( (10 ** 2, 10 ** 2) ), 0

	for key, (divergence, identity) in dict_div_id.items():
		
		divergence 	= int(round(divergence, 0))
		identity 	= int(round(identity, 0))
		
		div_ind = int(divergence)
		id_ind 	= int(identity)
		
		try:
			my_array[id_ind, div_ind] += 1
			total += 1
		except:
			pass;			# no identity information
				
	return my_array	
"""	
	
		
# SoDA alignment

def get_soda_summary(s):
	s 		= s.strip().split(sep)[ : -1]
	sid 	= s[0][s[0].rindex(" ") + 1 :]
	score 	= s[1].split("=")[1].strip()

	return sid, score


def get_soda_sequence(s):
	s = s[ 22 : ].rstrip()		# remove the description at the beginning
	return "".join([x for x in s if x!= " "]).upper()



def get_soda_desc(s):
	s = s.strip().split("=")

	func = s[0].split(" ")[1]
	
	
	if func == "Functional":
		reason, is_func = "", True
	else:
		reason, is_func = s[1], False
		
	return is_func, reason
	


def get_soda_seq(f):
	# return seq_id, aa seq, score
	
	handle 					= open(f, "rU")
	sid, score 				= "", ""
	is_func, reason			= "", ""
	nucleotides, protein 	= "", ""
	put_germ, germ			= "", ""
	
	for line in handle:
		#line = line.rstrip()
		
		if line.find("Sequence name") >= 0:
			sid, score = get_soda_summary(line)
		
		elif line.find("Functionality") >= 0:
			is_func, reason = get_soda_desc(line)
			
		elif line.find("protein") >= 0:
			protein += get_soda_sequence(line)
			
		elif line.find("input seq.") >= 0:
			nucleotides += get_soda_sequence(line)
			
		elif line.find("put. germline") >= 0:
			put_germ += get_soda_sequence(line)
			
		elif line.find("V ") >= 0:
			germ += get_soda_sequence(line)
	
	if nucleotides.find("N") >= 0:
		nucleotides = replace_N(put_germ, germ, nucleotides)

	if nucleotides.find("-") >= 0:
		nucleotides = remove_dash(nucleotides)

	nucleotides, protein = paired_stop_codon_trim(nucleotides, protein)
	
		
			
	if nucleotides.find("N") >= 0:
		print "N inside: %s" %sid, nucleotides.index("N"), len(nucleotides)
		print nucleotides
		print
	
	if nucleotides.find("-") >= 0:
		print "- inside: %s" %sid
		print
	
	protein_len = len(protein)
	if 0 <= protein.find("*") < protein_len - 1:
		print "STOP codon inside: %s; length: %d; stop@  %d" %(sid, protein_len, protein.index("*"))
		print
	
	return sid, score, is_func, reason, protein, nucleotides


def paired_stop_codon_trim(nuc_seq, aa_seq):
	
	aa_len, nuc_len = len(aa_seq), len(nuc_seq)
	
	if aa_len * 3 - nuc_len != 0:
		aa_len -= 1
		nuc_len = aa_len * 3
		
	nuc_seq, aa_seq = nuc_seq[ : nuc_len], aa_seq[ : aa_len]
	
	if aa_seq.find("*") < 0:
		return nuc_seq, aa_seq
		
	ind = aa_seq.index("*")
	
	aa_seq = aa_seq[ : ind]
	nuc_seq = nuc_seq[ : 3 * ind]
	
	
	return nuc_seq, aa_seq

	

def load_read_assigned(align_file, germ, min_len):
	print "generating read germ alignment from file: %s; germ: %s; min len: %d" %(align_file, germ, min_len)
	reader, total, dict_read_germ = csv.reader(open(align_file, "rU"), delimiter = sep), 0, dict()
	reader.next()
	for row in reader:
		qid, sid, strand, sstart, send = row[0].strip(), row[1].strip(), row[-3], int(row[8]), int(row[9])
		aln_len = abs(sstart - send + 1)
		
		if sid == germ and aln_len >= min_len:
			dict_read_germ[qid] = strand
	
	print "%d in total" %len(dict_read_germ)
	return dict_read_germ

"""
def load_read2germ_assignment(align_file , germ):
	print "loading reads assigned to germline: %s from file: %s..." %(germ, align_file)
	result, total 	= dict(), 0
	reader			= csv.reader(open(align_file, "rU"), delimiter = sep)
	for row in reader:
		qid, sid = row[0].strip(), row[1].strip()
		if sid == germ:
			result[qid] = [nan] * 2
			total += 1
	
	print "%d reads loaded..." %total
	return result
"""
	

def load_read_assigned2germ(align_file, germs):
	print "loading reads assigned to germline from file: %s..." %(align_file)
	result, total 	= set(), 0
	reader 			= csv.reader(open(align_file, "rU"), delimiter = sep)
	for row in reader:
		qid, sid = row[0].strip(), row[1].strip()
		if sid in germs:
			result.add(qid)
			total += 1

	print "%d reads loaded..." %total
	return result
		

def generate_read_native(folder, native, max_start, min_end):
	"""
	generate reads that were aligned to native with max start and min end (refering to QVQ/VSS)
	"""
	
	infiles = glob.glob("%s/*.txt" %folder)	
	for infile in infiles:
		print "file: %s... native: %s; max_start: %d; min_end: %d" %(infile, native, max_start, min_end)
		reader = csv.reader(open(infile, "rU"), delimiter = sep)
		for row in reader:
			qid, sid = row[0], row[1]
			if sid == native:
				sstart, send = int(row[-4]), int(row[-3])
				sstart, send = min(sstart, send), max(sstart, send)
				if sstart <= max_start and send >= min_end:
					yield qid
		#break

def calculate_levenshtein_distance(s1, s2):
	l1, l2 = len(s1) + 1, len(s2) + 1
	dist_array = zeros((l1, l2))
	for i in range(l1):
		dist_array[(i, 0)] = i
	for i in range(l2):
		dist_array[(0, i)] = i
	#print dist_array

	for i in range(1, l1):
		for j in range(1, l2):
			if s1[i - 1] == s2[j - 1]:
				dist_array[(i, j)] = dist_array[(i-1, j-1)]
			else:
				dist_array[(i, j)] = min(dist_array[(i-1, j)] + 1, dist_array[(i, j-1)] + 1, dist_array[(i-1, j-1)] + 1)

	#print dist_array
	return dist_array[l1 - 1, l2 - 1]

# SoDA2.0 alinment

class SoDA2:
	""""""
	def __init__(self):
		self.gene_id 	= ""
		self.is_func	= False
		self.v			= ""
		self.d			= ""
		self.j			= ""
		self.score		= ""
		
		self.input_seq	= ""
		self.put_germ	= ""
		self.key		= ""
		self.input_aa	= ""
		self.germ_aa	= ""
		

def get_soda_seq_id(s):
	""" Sequence Name: 124822_3 """
	return s[s.rindex(" ") + 1 :].strip()
	

def get_soda_func_reason(s):
	""" Human Ig heavy sequence                 
	Functionality : Non-Functional                  
	Reason : Stop Codon in Input Sequence 
	"""
	
	s = s.strip().split(sep)
	is_func, reason = True, ""
	
	if s[3].find("Non-Functional") > 0:
		is_func = False
	
	if not is_func:
		reason = s[6][s[6].index(":") + 2 :]
	
	return is_func, reason


def get_vdj(line):
	""" V Gene: 1~2*02  D Gene: 3-3*02  J Gene: 1*01 """
	def retrieve_value(s):
		ind = s.rindex(" ") + 1
		return s[ind : ]
		
	return tuple(map(retrieve_value, line.split(sep)))
	
def parse_soda2_seq(line):
	""" Input Seq   C A G G C A C A C C T G G A G C A """
	
	return "".join(line[12 :].split(" "))

def parse_soda2_aa(line):
	codon_sep = "     "
	return "".join(line[12 :].strip().split(codon_sep))


def vet_seq(my_soda2, f, sid):
	""" curate the sequence """
	
	def remove_insertion(input_codon, germ_codon):
		z = zip(input_codon, germ_codon)
		return "".join([x for x, y in z if y != "~"])
		
	def replace_nucleotide(input_codon, germ_codon, nucs):
		z, result = zip(input_codon, germ_codon), []
		for x, y in z:
			if x in nucs:
				x = y
			result.append(x)
		return "".join(result)
		
	def alter_stop(org, ref):
		z, result = zip(org, ref), []
		for (x, y) in z:
			if x != y and y != "-" and y != "~":
				x = y
			result.append(x)
			
		return "".join(result)
		
	input_seq 	= "".join(my_soda2.input_seq)
	put_germ	= "".join(my_soda2.put_germ)
	keys		= "".join(my_soda2.key)
	input_aa	= "".join(my_soda2.input_aa)
	germ_aa		= "".join(my_soda2.germ_aa)
	
	input_aa, germ_aa = get_soda2_codon_list(input_aa), get_soda2_codon_list(germ_aa)
	aa_len = len(input_aa)
	
	input_seq, put_germ = remove_odd(input_seq)[ : 3 * aa_len], remove_odd(put_germ)[ : 3 * aa_len]
	keys				= remove_odd(keys)[ : 3 * aa_len]
	
	input_seq = replace_nucleotide(input_seq, put_germ, "N-")
	
	seq_list, germ_list = [], []
	
	result, start = [], 0
	for ind, aa in enumerate(input_aa):
		codon_len 	= len(aa)
		seq_len 	= codon_len + 2
		end = start + seq_len
		
		input_codon = input_seq[start : end]
		germ_codon	= put_germ[start : end]
		
		seq_list.append(input_codon)
		germ_list.append(germ_codon)
		
		start = end			# record end position
		
		if codon_len == 1:
			if input_codon in stop_codons:
				new_codon = alter_stop(input_codon, germ_codon)
			else:
				new_codon = input_codon
				
		elif codon_len == 2:
				new_codon = remove_insertion(input_codon, germ_codon)
				
				if len(new_codon) != 3:	# > 1 insertion, take first 3 input
					new_codon = input_codon[ : 3]
				

		elif codon_len > 2:			# more than one insertion
			new_codon = remove_insertion(input_codon, germ_codon)
			if len(new_codon) > 3:
				new_codon = new_codon[ : 3]
			elif len(new_codon) < 3:
				new_codon = input_codon[ : 3]
		else:
			print "^&^   ^&^%s^&^   ^&^" %aa
			new_codon = ""
		
		if len(new_codon) == 3:
			result.append(new_codon)
		
		if new_codon.find("N") >= 0 or new_codon.find("-") >= 0 or new_codon.find("~") >= 0:
			break
		
	result = "".join(result)
	return result
	
	

	
def remove_odd(s):			# odd coz python start from "0"
	return s[ :: 2]	


def get_soda2_codon_list(s):
	return remove_odd(s).strip().split("  ")
	

def parse_SoDA2_alignments(folder):
	infiles 	= get_files_format(folder, "txt")
	dict_soda 	= dict()

	my_soda2 = SoDA2()
	for infile in infiles:		
		reader = open(infile, "rU")
		for line in reader:
			if line.find("Sequence Name") >= 0:		
				
				if my_soda2.gene_id != "":		# old object, need to store in dict
					infile, my_soda2.gene_id
					dict_soda[my_soda2.gene_id] = (vet_seq(my_soda2, infile, my_soda2.gene_id), infile)
					
				my_soda2 = SoDA2()		# create a new object
				my_soda2.gene_id = get_soda_seq_id(line)
				
			
			if line.find("Functionality") >= 0:
				my_soda2.isfunc, my_soda2.reason = get_soda_func_reason(line.strip())
				
			if line.find("V Gene") >= 0:
				my_soda2.v, my_soda2.d, my_soda2.j = get_vdj(line.strip())
			
			if line.startswith("Input Seq"):
				my_soda2.input_seq = my_soda2.input_seq + line[12 : -1]
			
			if line.startswith("Put Germ"):
				my_soda2.put_germ = my_soda2.put_germ + line[12 : -1]
				
			if line.startswith("Key"):
				my_soda2.key = my_soda2.key + line[12 : -1]
				
			if line.startswith("InputAA"):
				my_soda2.input_aa = my_soda2.input_aa + line[12 : -1]
				
			if line.startswith("Germ AA"):
				my_soda2.germ_aa = my_soda2.germ_aa + line[12 : -1]
				
	dict_soda[my_soda2.gene_id] = (vet_seq(my_soda2, infile, my_soda2.gene_id), infile)		# last entry		
	
	return dict_soda
				
			
#
# -- END -- alignment file methods 
#

def get_germ_fam(s):
	return s[ : s.index("*")]


#
# -- BEGIN -- numbers, list, numpy
#

def my_float_range_inclusive(start, end, step):
	# return a list with given step, both end included
	result = [start]
	
	while (end - start) > 0.001:
		start += step
		result.append(round(start, 2))
		
	return result

def is_any_identity_higher(l, min_id, indice):
	for index in indice:
		if l[index] >= min_id:
			return True
	
	return False

def is_higher(l, min_id):
	return sum([x > min_id for x in l])

#
# -- END -- numbers, list, numpy
#





#
# -- BEGIN -- Methods from pevious mytools
#


def swapExt(s):
	if s == ".txt":
		return ".tab"
	if s == ".tab":
		return ".txt"
		

def get_title_content(f_csv):
	reader, master_list = csv.reader(open(f_csv, "rU"), delimiter = sep), []

	title = reader.next()
	for i in range(len(title)):
		master_list.append([])

	for ind, row in enumerate(reader):
		
		row = map(float, row)
		for ind, item in enumerate(row):
			master_list[ind].append(item)
	return title, master_list
	


def binList(l, bin_size):
	
	result = []
	if len(l) % bin_size <> 0:
		raise Exception("Please select bin size appropriately")
		
		
	if not validateList(l):
		raise Exception("Invalid value in the list: None")
		
		
	for i in range(0, len(l) / bin_size):
		result.append(sum(l[i * bin_size : (i + 1) * bin_size]))
		
		
	return result



def smoothList(l, smooth_factor):
	#smooth the list value;
	#calculate the sum of surrounding smooth_factor values together with ith value

	result, l_len = [], len(l)

	if not validateList(l):
		raise Exception("invalid value in the list: None")

	for i in range(0, l_len):
		l_lim = max(0, i - smooth_factor/2)
		r_lim = min(l_len, i + smooth_factor/2 + 1)
		result.append( float(sum(l[l_lim : r_lim])) / float(r_lim - l_lim))

	return result



def validateList(l):
	# validate the list according to its type

	return not (None in set(l))



def load_divid(native, is_light):
	infile = ""
	if is_light == 0:
		infile = "%s/db/native/Heavy_divid_%s.txt" %(HOME_FOLDER, native)   #"Heavy_divid_%s.txt" %native
		xs, ys, natives = [], [], []
		reader = csv.reader(open(infile, "rU"), delimiter = sep)
		reader.next()
		for row in reader:
			natives.append(row[0])
			xs.append(float(row[1]))
			ys.append(float(row[2]))
		#print natives
		return natives, xs, ys
	
#
# -- END -- Methods from pevious mytools
#

def load_phylo_native_CDRH3():
	nat_cdr3_file, result = "%s/db/native/vrc_IgHV_CDRH3.fa" %HOME_FOLDER, dict()
	for Ab in SeqIO.parse(open(nat_cdr3_file, "rU"), "fasta"):
		if Ab.id in PHYLO_NATIVES:
			result[Ab.id] = Ab
	print "%d native CDRH3 loaded..." %len(result)
	
	#print set(PHYLO_NATIVES).difference(set(result.keys()))
	
	return result

def load_all_native_CDRH3():
	nat_cdr3_file, result = "%s/db/native/vrc_IgHV_CDRH3.fa" %HOME_FOLDER, dict()
	for Ab in SeqIO.parse(open(nat_cdr3_file, "rU"), "fasta"):
		if len(Ab.seq) > 0:
			result[Ab.id] = Ab
	print "%d native CDRH3 loaded..." %len(result)
	
	#print set(PHYLO_NATIVES).difference(set(result.keys()))
	
	return result
	
	
