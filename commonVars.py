#!/usr/bin/env python3

"""
commonVars.py

Created by Zhenhai Zhang on 2011-04-06.
Renamed by Chaim A Schramm on 2015-02-10.

Copyright (c) 2011-2019 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
"""

import sys, os, csv, shutil, re, glob, string, time, random
try:
	from SONAR.paths import *
except ImportError:
	sys.exit("Can't find paths.py. Have you run setup yet?")
from Bio.Data import CodonTable

sep = "\t"
linesep = os.linesep


# patterns
pat_nuc_cxrk = re.compile("TG[T|C]...(CG.|AG[A|G]|AA[A|G])", flags=re.I)


# programs (many are specified in paths.py after install)
igphyml        = "%s/third-party/igphyml_blas_omp"       % SCRIPT_FOLDER
igphyml_slow   = "%s/third-party/igphyml_no_libraries"   % SCRIPT_FOLDER
reconstruct    = "%s/third-party/ancReconstructHLP17.pl" % SCRIPT_FOLDER


#control variables
SUPPORTED_SPECIES = { "human":"HU", "rhesus":"RH" }
LOCUS_LIST        = [ 'H', "K", "L", "KL", "HKL" ]

# databases
HU_VH_DB	= "%s/germDB/IgHV.fa"   %SCRIPT_FOLDER
HU_VK_DB	= "%s/germDB/IgKV.fa"   %SCRIPT_FOLDER
HU_VL_DB	= "%s/germDB/IgLV.fa"   %SCRIPT_FOLDER
HU_VKL_DB	= "%s/germDB/IgKLV.fa"  %SCRIPT_FOLDER
HU_VHKL_DB	= "%s/germDB/IgHKLV.fa" %SCRIPT_FOLDER

HU_JH_DB	= "%s/germDB/IgHJ.fa"   %SCRIPT_FOLDER
HU_JK_DB	= "%s/germDB/IgKJ.fa"   %SCRIPT_FOLDER
HU_JL_DB	= "%s/germDB/IgLJ.fa"   %SCRIPT_FOLDER
HU_JKL_DB	= "%s/germDB/IgKLJ.fa"  %SCRIPT_FOLDER
HU_JHKL_DB	= "%s/germDB/IgHKLJ.fa" %SCRIPT_FOLDER

HU_DH_DB   = "%s/germDB/IgHD.fa"  %SCRIPT_FOLDER

HU_CH_DB	= "%s/germDB/IgHC_CH1.fa" %SCRIPT_FOLDER
HU_CK_DB	= "%s/germDB/IgKC.fa" %SCRIPT_FOLDER
HU_CL_DB	= "%s/germDB/IgLC.fa" %SCRIPT_FOLDER
HU_CKL_DB	= "%s/germDB/IgKLC.fa" %SCRIPT_FOLDER
HU_CHKL_DB	= "%s/germDB/IgHKLC.fa" %SCRIPT_FOLDER


RH_VH_DB	= "%s/germDB/IgHV_BU_DD.fasta"   %SCRIPT_FOLDER
RH_VK_DB	= "%s/germDB/IgKV_BU_DD.fasta"   %SCRIPT_FOLDER
RH_VL_DB	= "%s/germDB/IgLV_BU_DD.fasta"   %SCRIPT_FOLDER
RH_VKL_DB	= "%s/germDB/IgKLV_BU_DD.fasta"  %SCRIPT_FOLDER
RH_VHKL_DB	= "%s/germDB/IgHKLV_BU_DD.fasta" %SCRIPT_FOLDER

RH_JH_DB	= "%s/germDB/IgHJ_BU_DD.fasta"   %SCRIPT_FOLDER
RH_JK_DB	= "%s/germDB/IgKJ_BU_DD.fasta"   %SCRIPT_FOLDER
RH_JL_DB	= "%s/germDB/IgLJ_BU_DD.fasta"   %SCRIPT_FOLDER
RH_JKL_DB	= "%s/germDB/IgKLJ_BU_DD.fasta"  %SCRIPT_FOLDER
RJH_HKL_DB	= "%s/germDB/IgHKLJ_BU_DD.fasta" %SCRIPT_FOLDER

RH_DH_DB   = "%s/germDB/IgHD_BU_DD.fasta"  %SCRIPT_FOLDER

RH_CH_DB	= "%s/germDB/IgHC_CH1_BU_DD.fasta" %SCRIPT_FOLDER
RH_CK_DB	= "%s/germDB/IgKC_BU_DD.fasta" %SCRIPT_FOLDER
RH_CL_DB	= "%s/germDB/IgLC_BU_DD.fasta" %SCRIPT_FOLDER
RH_CKL_DB	= "%s/germDB/IgKLC_BU_DD.fasta" %SCRIPT_FOLDER
RH_CHKL_DB	= "%s/germDB/IgHKLC_BU_DD.fasta" %SCRIPT_FOLDER


#J motifs
HU_JH_MOTIF   = "TGGGG"
HU_JK_MOTIF   = "TT[C|T][G|A]G"
HU_JL_MOTIF   = "TT[C|T][G|A]G"
HU_JKL_MOTIF  = "TT[C|T][G|A]G"
HU_JHKL_MOTIF = "(TGGGG|TT[C|T][G|A]G)"

RH_JH_MOTIF   = "TGGGG"
RH_JK_MOTIF   = "TT[C|T][G|A]G"
RH_JL_MOTIF   = "TT[C|T][G|A]G"
RH_JKL_MOTIF  = "TT[C|T][G|A]G"
RH_JHKL_MOTIF = "(TGGGG|TTCGG|TT..GAA|TTTGGC|TTCTGT)"


CMD_BLAST           = "%s %s %s -query %s -out %s -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand' -gapopen 5 -gapextend 2 -penalty -1 -reward 1 -evalue 1e-3 -max_target_seqs 10 -word_size %d"
V_BLAST_WORD_SIZE   = 7
J_BLAST_WORD_SIZE   = 5
BLAST_OTHER_OPTIONS = "-evalue 1e-3"

PBS_STRING = "\
#!/bin/bash\n\
#$ -N %s		# job name\n\
#$ -l mem=%s,time=%s   # resource requests\n\
#$ -cwd				# use current directory as job status output\n\
#$ -o /dev/null			# use sane outputs for array jobs\n\
#$ -e /dev/null\n\
\n\
%s\n"


CMD_BLASTCLUST	= "/ifs/home/c2b2/bh_lab/shares/blast/current/ia32-linux/bin/blastclust -p F -L .9 -S 95 -i %s -o %s"	#pF: nucleotide; L.9: 90%[coverage]  S: Identities


PARSED_BLAST_HEADER = ["qid", "sid", "identity", "align_len", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "score", "strand","other_sids"]
PARSED_BLAST_HEADER_VERBOSE = ["query_id", "sbjct_id", "strand", "evalue", "score", "identities", "gaps", "aln_len",
								"query_start", "query_end", "query_len", "sbjct_start", "sbjct_end", "aln_query", "aln_sbjct"]


#make a codon table that can handle gaps
#start by getting the standard codon table
table = CodonTable.standard_dna_table.forward_table
#add gaps
for c1 in ["A", "C", "G", "T", "N"]:
	table["%s--"%c1] = "X"
	table["-%s-"%c1] = "X"
	table["--%s"%c1] = "X"
	for c2 in ["A", "C", "G", "T", "N"]:
		table["%s%s-"%(c1,c2)] = "X"
		table["-%s%s"%(c1,c2)] = "X"
		table["%s-%s"%(c1,c2)] = "X"
table["---"]="-"
#now register is and export
CodonTable.register_ncbi_table(name='gapped',alt_name="CAS0",id=99,table=table, stop_codons=['TAA', 'TAG', 'TGA', ], start_codons=['TTG', 'CTG', 'ATG', ] )
GAPPED_CODON_TABLE=CodonTable.ambiguous_dna_by_name["gapped"]
