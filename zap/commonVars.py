#!/usr/bin/env python
# encoding: utf-8
"""
common_info.py

Created by Zhenhai Zhang on 2011-04-06.
Copyright (c) 2011 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
"""

import sys, os, csv, shutil, re, glob, pickle, string, time, random, commands

sep = "\t"
linesep = os.linesep


# patterns
# heavy chain
pat_nuc_cxrk = re.compile("TG[T|C]...(CG.|AG[A|G]|AA[A|G])")
pat_nuc_wgxg = re.compile("TGGGG....GG.")
pat_nuc_ast  = re.compile("GC.(TC.|AG[T|C])AC.")
pat_nuc_vss  = re.compile("GT.(TC.|AG[A|C])(TC.|AG[A|C])")

pat_pro_cxrk = re.compile("C.[R|K]")
pat_pro_wgxg = re.compile("WG.G")
pat_pro_ast	 = re.compile("AST")

# light chain
pat_nuc_yxc  = re.compile("TA[T|C]...TG[T|C]")
pat_nuc_cqq  = re.compile("TG[T|C]CA[A|G]CA[A|G]")
pat_nuc_fgxg = re.compile("TT[T|C]CG....CG.")
pat_nuc_gxg  = re.compile("GG....GG.")

v_pattern = re.compile("V_gene=(IG.*?)[,| ]")
#v_pattern = re.compile("(IG.*)")



#
# ===START=== folders, databases
#

HOME_FOLDER    = "/ifs/scratch/c2b2/ls_nih/cs3037"
clustalw       = "/ifs/home/c2b2/bh_lab/shares/clustalw/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2"
clustal	       = "/ifs/home/c2b2/bh_lab/cs3037/bin/clustalo-1.1.0-linux-64"
blastclust_cmd = "/ifs/home/c2b2/bh_lab/shares/blast/current/ia32-linux/bin/blastclust"
blastall_cmd   = "/ifs/home/c2b2/bh_lab/shares/blast/current/ia32-linux/bin/blastall"


#assert os.path.isfile(clustal), "Clustal W missing"

# databases	
VH_DB 		= "%s/db/germline/IgHV.fa" %HOME_FOLDER
VK_DB		= "%s/db/germline/IgKV.fa" %HOME_FOLDER
VL_DB		= "%s/db/germline/IgLV.fa" %HOME_FOLDER
VKL_DB		= "%s/db/germline/IgKLV.fa" %HOME_FOLDER

JH_DB 			= "%s/db/germline/IgHJ.fa" %HOME_FOLDER
JK_DB 			= "%s/db/germline/IgKJ.fa" %HOME_FOLDER
JL_DB 			= "%s/db/germline/IgLJ.fa" %HOME_FOLDER
JKL_DB 			= "%s/db/germline/IgKLJ.fa" %HOME_FOLDER

DH_DB                   = "%s/db/germline/IgHD.fa" %HOME_FOLDER
CH_DB                   = "%s/db/germline/IgHC.fa" %HOME_FOLDER


dict_vgerm_db = {
	0: VH_DB,
	1: VK_DB,
	2: VL_DB,
	3: VKL_DB
}

dict_jgerm_db = {
	0: JH_DB,
	1: JK_DB,
	2: JL_DB,
	3: JKL_DB,
}

# amino acid databases
GERM_HEAVY_AA	= "%s/db/germline/IgHV_protein.fasta" %HOME_FOLDER
NAT_HEAVY_AA	= "%s/db/native/vrc_IgHV_protein.fasta" %HOME_FOLDER


ALL_FOLDERS = ["work/1-blast", "work/2-clustal", "work/3-phylogeny", "work/1-blast/vgene", "work/1-blast/jgene", "output/sequences", "output/tables", "output/plots", "output/logs"]


#
# ===END=== folders, databases
#


#
# ===START=== PBS
#

CMD_BLASTALL            = "%s -p blastn -m 8 %s -d %s -i %s -o %s" % (blastall_cmd,"%s","%s","%s","%s")
BLAST_V_OPTIONS 	= " -J T -G 5 -E 2 -q -1 -r 1 -W 7 -b 5 -v 1 -e 1e-10 "
BLAST_J_OPTIONS 	= " -J T -G 5 -E 2 -q -1 -r 1 -W 5 -b 5 -v 1 -e 1e-3 "
BLAST_OTHER_OPTIONS = "-e 1e-3"

PBS_STRING = "\
#!/bin/bash\n\
#$ -t 1-%d			# array task\n\
#$ -N %s		# job name\n\
#$ -l mem=%s,time=%s	# resource requests\n\
#$ -cwd				# use current directory as job status output\n\
#$ -o /dev/null			# use sane outputs for array jobs\n\
#$ -e /dev/null\n\
\n\
%s\n"


CMD_BLASTCLUST	= "/ifs/home/c2b2/bh_lab/shares/blast/current/ia32-linux/bin/blastclust -p F -L .9 -S 95 -i %s -o %s"	#pF: nucleotide; L.9: 90%[coverage]  S: Identities 


#
# ===END=== PBS
#


#
# ===START=== file headers
#

PARSED_BLAST_HEADER = ["qid", "sid", "identity", "align_len", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "score", "strand", "slen","other_sids"]
PARSED_BLAST_HEADER_VERBOSE = ["query_id", "sbjct_id", "strand", "evalue", "score", "identities", "gaps", "aln_len", 
								"query_start", "query_end", "query_len", 
								"sbjct_start", "sbjct_end", "sbjct_len", 
								"aln_query", "aln_sbjct"]



#
# ===END=== file headers
#


DICT_H_NATIVES = {
	"IGHV1-2*02" : ("VRC01", "VRC03", "VRC04")
}

LIST_NATIVES_BAK   = ["VRC01", "VRC03", "VRC06", "VRC06B", "VRC07", "VRC08", "VRC02"]

#LIST_NATIVES  = ['VRC01', 'VRC01B', 'VRC02', 'VRC03', 'VRC03B', 'VRC03C', 'VRC06', 'VRC06B', 'VRC07', 'VRC07B', 'VRC07C', 'NIH45-46', 'VRC08', 
#				'VRC08B', 'VRC-PG04', 'VRC-PG04B', 'VRC17', 'VRC18', 'VRC-PG19', 'VRC-PG19B', 'VRC-PG20', 'VRC-PG20B', 'VRC-CH30', 'VRC-CH31', 
#				'VRC-CH32', 'VRC-CH33', 'VRC-CH34', '3BNC60', '3BNC117', '12A12', '12A21', "CH01", "CH02", "CH03", "CH04", "CH103", "CH104", 
#				"CH105", "CH106"]

LIST_NATIVES   = ['VRC01','VRC03','VRC-PG04','VRC-PG04_cog','VRC-CH31','VRC-CH33','gVRC-d74_H3','gVRC-d74_H6','gVRC-d74_H12','gVRC-d74_H15']
LIST_VALID     = ['VRC01','VRC03','VRC06','VRC08','gVRC-d45-H3','gVRC-d45-H7','gVRC-d45-H15','gVRC-d45-H16','gVRC-d45-H17','gVRC-d45-H18','gVRC-d45-H19','gVRC-d45-H27','gVRC-d45-H29','gVRC-d45-H30','gVRC-d45-H31','gVRC-d45-H32','gVRC-d45-H33','gVRC-d45-H34']




LIST_L_NATIVES = ["VRC01", "VRC02", "VRC03", "VRC-PG04", "VRC-PG04B", "VRC06", "VRC06B", "PG16", "PG9", "VRC-CH30", "VRC-CH31", "VRC-CH32", 
				"CH01", "CH02", "CH03", "CH04", "CH103", "CH104", "CH105", "CH106", "VRC07B", "VRC07C", "VRC08", "12A21", "12A12", "12A21", 
				"NIH45-46", "1B2530", "1NC9", "8ANC131", "8ANC134"]

LIST_CHAVI_H = ["UCA","I1","I2","I3","I4","I7","I8","1AZCETI5","1A102RI6","1AH92U","CH103","CH104","CH106"]
LIST_CHAVI_L = ["UCA","53.I698R","53.II1CU","I1","I2","CH103","CH104","CH106"]
LIST_CH505 = ["CH103","CH104","CH106","CH186","CH187","CH188","CH200"]

LIST_CAP256_H = ["VRC26","VRC26b","VRC26c","VRC26d","VRC26e","VRC26f","VRC26g","VRC26h","VRC26i","VRC26j","VRC26k","VRC26m","UCA","I1","I2"]
LIST_CAP256_L = ["VRC26","VRC26b","VRC26c","VRC26d","VRC26e","VRC26f","VRC26g","VRC26h","VRC26i","VRC26j","VRC26k","VRC26m","UCA","I1","I2"]
LIST_CAPTEST = ["CAP256_Ab1","CAP256_Ab2","CAP256_Ab3","CAP256_Ab4"]

LIST_10E8_L = ["10E8"]
LIST_PGT141_H = ['PGT141','PGT142','PGT143','PGT144','PGT145']

LIST_NUSSEN	   = ['VRC01', 'VRC02', 'VRC03', 'VRC-PG04', 'VRC-PG04B', 'VRC06', 'VRC06B', "12A21", "12A12", "12A21", "NIH45-46", "1B2530", "1NC9", 
				"8ANC131", "8ANC134", "8ANC195"]
LIST_H_MONKEY  = ["GE121_HC", "GE125_HC", "GE136_HC", "GE137_HC", "GE140_HC", "GE143_HC", "GE147_HC", "GE148_HC"]

LIST_NIH57 = ["00076734","00130597","00082681","00182776"]

DICT_NATIVES = {
	0 : LIST_NATIVES,
	1 : LIST_L_NATIVES,
	2 : LIST_L_NATIVES, #for testing with lamba, too; original: LIST_H_MONKEY,
	3 : LIST_L_NATIVES, #for testing with lamba, too; original: []
	
	6 : LIST_NIH57,

	10 : LIST_CH505,
	11 : LIST_CH505,

	12 : LIST_CAP256_H,
	13 : LIST_CAP256_L,
	14 : LIST_CAPTEST,

	15 : LIST_10E8_L,
	16 : LIST_PGT141_H,
	18 : LIST_VALID,

	99: LIST_NUSSEN,
	100 : LIST_NUSSEN,
	1000: LIST_NATIVES
}

DICT_PHYLO_INFO = {
	30 	: (GERM_KCL_HEAVY, NAT_KCL_HEAVY, ["VH3.8", "VH3.18", "VH3.24"]), # (GERM FILE, NATIVE FILE, [GERM V GENE LIST])
	40 	: (GERM_KCL_HEAVY, NAT_KCL_HEAVY, ["VH4.11", "VH4.40", "VH4.57"]),
	0	: (GERM_HEAVY, NAT_HEAVY, ["IGHV1-2*02"])
}

DICT_PHYLO_GERM_NAT = {
	"VH3.8"		: ["GE121_HC", "GE125_HC"],
	"VH3.18"	: ["GE137_HC"],
	"VH3.24"	: ["GE147_HC"],
	"VH4.11"	: ["GE136_HC", "GE140_HC"],
	"VH4.40"	: ["GE143_HC"],
	"VH4.57"	: ["GE148_HC"],
	"IGHV1-2*02": ["VRC01", "VRC01B", "VRC02", "VRC03", "VRC03B","VRC03C","VRC06", "VRC06B", "VRC07", "VRC07B", "VRC07C","VRC08","VRC08B", 
					"NIH45-46", "VRC-PG04", "VRC-PG04B", "VRC-CH30", "VRC-CH31", "VRC-CH32", "VRC-CH33", "VRC-CH34", "3BNC60", 
					"3BNC117", "12A12", "12A21", "1NC9", "1B2530", "8ANC131"]		# no CDRH3 available for "8ANC134"
}



#PHYLO_NATIVES  = ["VRC01", "VRC02", "VRC03", "NIH45-46", "VRC-PG04", "VRC-PG04B", "VRC-CH30", "VRC-CH31", "VRC-CH32", "VRC-CH33", "VRC-CH34", "3BNC60", 
#				"3BNC117", "12A12", "12A21", "1NC9", "1B2530", "8ANC131", "8ANC134"]
PHYLO_NATIVES  = ["VRC01", "VRC02", "VRC03", "VRC06", "VRC07", "VRC08", "NIH45-46", "VRC-PG04", "VRC-CH30", "VRC-CH31",
		                "VRC-CH32", "VRC-CH33", "VRC-CH34", "3BNC60", "3BNC117", "12A12", "12A21"]

PHYLO_CDRH3 = ["NIH45-46","VRC17","CH02","VRC-CH30","VRC-PG20","VRC18","VRC26","VRC26j","VRC26i","10e8","PGT128","PGT145","PG9","8ANC195"]#,"CH103"]
PHYLO_CDRL3 = ["VRC03","NIH45-46","CH02","VRC-CH30","VRC26","VRC26j","VRC26i","10e8","PGT128","PGT145","PG9","8ANC195","CH103"]#,"CH31"]

PHYLO_CHAVI = ["UCA","I1","I2","I3","I4","I7","I8","1AZCETI5","1A102RI6","1AH92U","CH103","CH104","CH106"]
PHYLO_CH505 = ["CH103","CH104","CH106","CH186","CH187","CH188","CH200"]

PHYLO_CAP256 = ["VRC26","VRC26b","VRC26c","VRC26d","VRC26e","VRC26f","VRC26g","VRC26h","VRC26i","VRC26j","VRC26k","VRC26m"]

PHYLO_PGT141 = ['PGT141','PGT142','PGT143','PGT144','PGT145']

PHYLO_CINQUE1 = ['VRC01', 'NIH45-46', 'VRC03']
PHYLO_CINQUE2 = ['VRC13', 'VRC13f', 'VRC14']
PHYLO_CINQUE3 = ['VRC01', 'NIH45-46', 'VRC03', 'VRC13', 'VRC13f', 'VRC14']

DICT_PHYLO_NITIVES = {
	0 : PHYLO_NATIVES,
	1 : PHYLO_NATIVES, 
	2 : [],

	4 : PHYLO_CDRH3,
	5 : PHYLO_CDRL3,

	6 : PHYLO_NATIVES, #did NIH57 XD analysis with this set

	10 : PHYLO_CH505,
	11 : PHYLO_CH505,
	
	12 : PHYLO_CAP256,
	13 : PHYLO_CAP256,
	
	15: LIST_10E8_L,
	16: PHYLO_PGT141,

	99 : LIST_NUSSEN,
	100 : LIST_NUSSEN,

	101 : PHYLO_CINQUE1,
	102 : PHYLO_CINQUE2,
	103 : PHYLO_CINQUE3,

	1000: PHYLO_NATIVES

}

SET_PHYLO_NATIVES = set(PHYLO_NATIVES)
LEN_PHYLO_NATIVES = len(PHYLO_NATIVES)



nucleotides = ["A", "C", "G", "T"]

#
# ===START=== amino acid
#
MAX_PRIMER_MISMATCHES = 3
# Primers in VRC antibody sequencing - 
# ---VH1 gene amplification
# forward primers
DICT_HEAVY_FORWARD_PRIMERS = {
	"L_VH1"		: "ACAGGTGCCCACTCCCAGGTGCAG",
	"L_VH1_2"	: "GCAGCCACAGGTGCCCACTCC",
	"L_VH1_24"	: "CAGCAGCTACAGGCACCCACGC",
	"L_VH1_69"	: "GGCAGCAGCTACAGGTGTCCAGTCC"
}
# reverse primers
# the actual value is the reverse complement of the real primer
#
# IgM=mu: 		targets naive B cells
# IgG=gamma:	targets memory B Cells
DICT_HEAVY_REVERSE_PRIMERS = {
	"gamma"	: "ACCAAGGGCCCATCGGTC", #TTCCCC", 			# old "CCACCAAGGGCCCATCGGTCTTCCCCC",		# length: 27
			   
	"mu"	: "TCGTCTCCTGTGAGAA", #TTCCC"				# length: 21
	#"unknown"	: ""
}

DICT_HEAVY_REVERSE_PRIMERS_LEN = {
	"gamma" : 24, #27, 
	"mu"	: 21
}

# VK3 forward
DICT_LIGHT_FORWARD_PRIMERS = {
	"L_VK3"	: "CTCTTCCTCCTGCTACTCTGGCTCCCAG"
}
# masked coz we need to reverse complement
DICT_LIGHT_REVERSE_PRIMERS = {
	"CK494"	: "AGCAGGACAGCAAGGACAGCAC"
}

# Rhesus Macaque
DICT_RM_HEAVY_REVERSE_PRIMERS = {
	"gamma" : "CCACCAAGGGCCCATCGGTCTTCCCC",
	"mu"    : "TCGTCTCCTGTGAGAATTCCC"
}
DICT_RM_HEAVY_REVERSE_PRIMERS_LEN = {
	"gamma" : 26, 
	"mu"	: 21
}

#
# === END === amino acid
#




#
# ===START=== amino acid
#

# codons	
start_codons 	= (["ATG"])
stop_codons 	= (["TAA", "TAG", "TGA"])



dict_codon2aa = {

	"ATT" : ("I", "Isoleucine"),
	"ATA" : ("I", "Isoleucine"),
	"ATC" : ("I", "Isoleucine"),
	
	"CTT" : ("L", "Leucine"),
	"CTC" : ("L", "Leucine"),
	"CTA" : ("L", "Leucine"),
	"CTG" : ("L", "Leucine"),
	"TTA" : ("L", "Leucine"),
	"TTG" : ("L", "Leucine"),
	
	"GTT" : ("V", "Valine"),
	"GTC" : ("V", "Valine"),
	"GTA" : ("V", "Valine"),
	"GTG" : ("V", "Valine"),
	
	"TTT" : ("F", "Phenylalanine"),
	"TTC" : ("F", "Phenylalanine"),
	
	"ATG" : ("M", "Methionine"),
	
	"TGT" : ("C", "Cysteine"),
	"TGC" : ("C", "Cysteine"),
	
	"GCT" : ("A", "Alanine"),
	"GCC" : ("A", "Alanine"),
	"GCA" : ("A", "Alanine"),
	"GCG" : ("A", "Alanine"),
	
	"GGT" : ("G", "Glycine"),
	"GGC" : ("G", "Glycine"),
	"GGA" : ("G", "Glycine"),
	"GGG" : ("G", "Glycine"),
	
	"CCT" : ("P", "Proline"),
	"CCC" : ("P", "Proline"),
	"CCA" : ("P", "Proline"),
	"CCG" : ("P", "Proline"),
	
	"ACT" : ("T", "Threonine"),
	"ACC" : ("T", "Threonine"),
	"ACA" : ("T", "Threonine"),
	"ACG" : ("T", "Threonine"),
	
	"TCT" : ("S", "Serine"),
	"TCC" : ("S", "Serine"),
	"TCA" : ("S", "Serine"),
	"TCG" : ("S", "Serine"),
	"AGT" : ("S", "Serine"),
	"AGC" : ("S", "Serine"),
	
	"TAT" : ("Y", "Tyrosine"),
	"TAC" : ("Y", "Tyrosine"),
	
	"TGG" : ("W", "Tryptophan"),
	
	"CAA" : ("Q", "Glutamine"),
	"CAG" : ("Q", "Glutamine"),
	
	"AAT" : ("N", "Asparagine"), 
	"AAC" : ("N", "Asparagine"),
	
	"CAT" : ("H", "Histidine"), 
	"CAC" : ("H", "Histidine"),
	
	"GAA" : ("E", "Glutamic acid"),
	"GAG" : ("E", "Glutamic acid"),
	
	"GAT" : ("D", "Aspartic acid"),
	"GAC" : ("D", "Aspartic acid"),
	
	"AAA" : ("K", "Lysine"),
	"AAG" : ("K", "Lysine"),
	
	"CGT" : ("R", "Arginine"),
	"CGC" : ("R", "Arginine"),
	"CGA" : ("R", "Arginine"),
	"CGG" : ("R", "Arginine"),
	"AGA" : ("R", "Arginine"),
	"AGG" : ("R", "Arginine"),
	
	"TAA" : ("*", "Stop codons"),
	"TAG" : ("*", "Stop codons"),
	"TGA" : ("*", "Stop codons")

}



dict_aa2codon = {

	"I" : ("ATT", "ATC", "ATA"),
	"L" : ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
	"V" : ("GTT", "GTC", "GTA", "GTG"),
	"F" : ("TTT", "TTC"),
	"M" : ("ATG"),
	"C" : ("TGT", "TGC"),
	"A" : ("GCT", "GCC", "GCA", "GCG"),
	"G" : ("GGT", "GGC", "GGA", "GGG"),
	"P" : ("CCT", "CCC", "CCA", "CCG"),
	"T" : ("ACT", "ACC", "ACG", "ACA"),
	"S" : ("TCT", "TCC", "TCG", "TCA", "AGT", "AGC"),
	"Y" : ("TAT", "TAC"),
	"W" : ("TGG"),
	"Q" : ("CAA", "CAG"),
	"N" : ("AAT", "AAC"),
	"H" : ("CAT", "CAC"),
	"E" : ("GAA", "GAG"),
	"D" : ("GAT", "GAC"),
	"K" : ("AAA", "AAG"),
	"R" : ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
	"STOP" : ("TAA", "TAG", "TGA")

}



dict_aa2name = {

	"I" : "Isoleucine",
	"L" : "Leucine",
	"V" : "Valine",
	"F" : "Phenylalanine",
	"M" : "Methionine",
	"C" : "Cysteine",
	"A" : "Alanine",
	"G" : "Glycine",
	"P" : "Proline",
	"T" : "Threonine",
	"S" : "Serine",
	"Y" : "Tyrosine",
	"W" : "Tryptophan",
	"Q" : "Glutamine",
	"N" : "Asparagine",
	"H" : "Histidine",
	"E" : "Glutamic acid",
	"D" : "Aspartic acid",
	"K" : "Lysine",
	"R" : "Arginine",
	"STOP" : "Stop condons"
}

#
# ===END=== amino acid
#


#
# ===START=== IUB
#


IUB_CODE = {
	"A" : ["A"], 
	"C" : ["C"], 
	"G" : ["G"], 
	"T" : ["T"], 
	"R" : ["A", "G"], 
	"Y" : ["C", "T"], 
	"K" : ["G", "T"], 
	"M" : ["A", "C"], 
	"S" : ["G", "C"], 
	"W" : ["A", "T"], 
	"B" : ["C", "G", "T"], 
	"D" : ["A", "G", "T"], 
	"H" : ["A", "C", "T"], 
	"V" : ["A", "C", "G"], 
	"N" : ["A", "C", "G", "T"]
}


	
IUB_complement = {
	"A" : "T",
	"C" : "G", 
	"G" : "C", 
	"T" : "A",
	"R" : "Y", 
	"Y" : "R",
	"K" : "M", 
	"M" : "K",
	"S" : "W",
	"W" : "S",
	"B" : "V",
	"V" : "B",
	"D" : "H",
	"H" : "D",
	"N" : "N"
}


#
# ===END=== IUB
#
DICT_YR_CHAR = {
	2001 : "A",
	2002 : "B", 
	2003 : "C",
	2004 : "D",
	2005 : "E",
	2006 : "F",
	2007 : "G", 
	2008 : "H", 
	2009 : "I",
	2010 : "J",
	2011 : "K"
}

DICT_DATASET_HEAD = {
	"H010730G1": "A",
	"H010730H1": "a",
	"H020416G1": "B",
	"H060714G1": "C",
	"H060714H1": "c",
	"H070109G1": "D",
	"H070712G1": "E",
	"H080117G1": "F",
	"H080117H1": "f",
	"H080819G1": "G",
	"H080819H1": "g",
	"H090602G1": "H",
	"H091231G1": "I",
	"HR080117G1": "F"
}

DICT_DATASET_TITLE = {
	"H010730G1": "07/30/2001",
	"H010730H1": "07/30/2001",
	"H020416G1": "04/16/2002",
	"H060714G1": "07/14/2006",
	"H060714H1": "07/14/2006",
	"H070109G1": "01/09/2007",
	"H070712G1": "07/12/2007",
	"H080117G1": "01/17/2008",
	"H080117H1": "01/17/2008",
	"H080819G1": "08/19/2008",
	"H080819H1": "08/19/2008",
	"H090602G1": "06/02/2009",
	"H091231G1": "12/31/2009",
	"HR080117G1": "01/17/2008",
	"H950615G1" : "03/20/1995"
}


DICT_LIGHT_HEAD = {
	"L010730" : "A",
	"L020416" : "B", 
	"L060714" : "C", 
	"L070109" : "D", 
	"L070712" : "E", 
	"L080117" : "F", 
	"L080819" : "G", 
	"L090602" : "H", 
	"L091231" : "I",
	"L950615" : "O"

}

DICT_NAT_CDRH3 = {
	"NIH45-46" 	: 	(294, 342),
	"VRC17"  	:	(294, 324),
	"CH02"  	: 	(294, 366),
	"VRC-CH30"  	: 	(321, 360),
	"VRC-PG20"  	: 	(294, 333),
	"VRC18"  	: 	(294, 324),
	"VRC26"         :       (294, 405),
	"VRC26j"        :       (263, 368),
	"VRC26i"        :       (291, 396),
	"CH103"         :       (291, 330),
	"10e8"  	: 	(300, 360),
	"PGT128"        :       (315, 372),
	"PGT145"        :       (294, 387),
	"PG9"           :       (291, 375),
	"8ANC195"       :       (303, 363),
}


DICT_NAT_CDRL3 = {
	"VRC03"		:	(261, 276),
	"VRC-CH30"	:	(264, 279),
	"CH31"      	:	(264, 279),
	"CH02"		:	(267, 294),
	"NIH45-46"	:	(258, 273),
	"CH103"         :       (252, 282),
	"VRC26"         :       (267, 303),
	"VRC26i"        :       (267, 291),
	"VRC26j"        :       (267, 291),
	"PGT128"        :       (255, 285),
	"10e8"          :       (261, 297),
	"PGT145"        :       (279, 306),
	"PG9"           :       (270, 300),
	"8ANC195"       :       (267, 294),
}

'''
DICT_NAT_CDRL3 = {
	"VRC01"		:	(258, 273),
	"VRC02"		:	(258, 273),
	"VRC03"		:	(261, 276),
	"VRC-PG04"	:	(258, 273),
	"VRC-PG04B"	:	(258, 273),
	"VRC-CH30"	:	(264, 279),
	"VRC-CH31"	:	(264, 279),
	"VRC-CH32"	:	(264, 279),
	"CH01"		:	(267, 282),
	"CH02"		:	(267, 282),
	"CH03"		:	(267, 282),
	"CH04"		:	(267, 282),
	"VRC06"		:	(261, 276),
	"VRC06B"	:	(261, 276),
	"VRC07B"	:	(258, 273),
	"VRC07C"	:	(258, 273),
	"NIH45-46"	:	(258, 273),
	"VRC08"		:	(267, 282)	
}
'''

DICT_RM_CDRH3 = {
	"GE121_HC"	: (285, 351),
	"GE125_HC"	: (285, 351),
	"GE136_HC"	: (285, 354),
	"GE137_HC"	: (291, 369),
	"GE140_HC"	: (285, 354),
	"GE143_HC"	: (288, 366),
	"GE147_HC"	: (285, 345),
	"GE148_HC"	: (288, 357)
}

DICT_NAT_CDRH1 = {
        "UCA"     :       (90, 105),
        "I1"      :       (90, 105),
        "I2"      :       (90, 105),
        "I3"      :       (90, 105),
        "I4"      :       (90, 105),
        "I7"      :       (90, 105),
        "I8"      :       (90, 105),
        "1AZCETI5":       (90, 105),
        "1AH92U"  :       (90, 105),
        "1A102RI6":       (90, 105),
        "CH103"         :       (90, 105),
        "CH104"         :       (90, 105),
        "CH106"         :       (90, 105)
}

DICT_NAT_CDRL1 = {
        "UCA"     :       (66,99),
        "I1"      :       (66,99),
        "I2"      :       (66,99),
        "CH103"         :       (66,90),
        "CH104"         :       (66,99),
        "CH106"         :       (66,99)
}


DICT_NAT_CDRH2 = {
	"UCA"     :       (147, 195),
	"I1"      :       (147, 195),
	"I2"      :       (147, 195),
	"I3"      :       (147, 195),
	"I4"      :       (147, 195),
	"I7"      :       (147, 195),
	"I8"      :       (147, 195),
	"1AZCETI5":       (147, 195),
	"1AH92U"  :       (147, 195),
	"1A102RI6":       (147, 195),
	"CH103"         : (147, 195),
	"CH104"         : (147, 195),
	"CH106"         : (147, 195)
}
