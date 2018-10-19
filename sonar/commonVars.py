#!/usr/bin/env python3
# encoding: utf-8
"""
commonVars.py

Created by Zhenhai Zhang on 2011-04-06.
Renamed by Chaim A Schramm on 2015-02-10.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
"""

import sys, os, csv, shutil, re, glob, string, time, random
from sonar.paths import *

sep = "\t"
linesep = os.linesep


# patterns
pat_nuc_cxrk = re.compile("TG[T|C]...(CG.|AG[A|G]|AA[A|G])", flags=re.I)


# databases	
VH_DB	= "%s/germDB/IgHV.fa"   %SCRIPT_FOLDER
VK_DB	= "%s/germDB/IgKV.fa"   %SCRIPT_FOLDER
VL_DB	= "%s/germDB/IgLV.fa"   %SCRIPT_FOLDER
VKL_DB	= "%s/germDB/IgKLV.fa"  %SCRIPT_FOLDER
VHKL_DB	= "%s/germDB/IgHKLV.fa" %SCRIPT_FOLDER

JH_DB	= "%s/germDB/IgHJ.fa"   %SCRIPT_FOLDER
JK_DB	= "%s/germDB/IgKJ.fa"   %SCRIPT_FOLDER
JL_DB	= "%s/germDB/IgLJ.fa"   %SCRIPT_FOLDER
JKL_DB	= "%s/germDB/IgKLJ.fa"  %SCRIPT_FOLDER
JHKL_DB	= "%s/germDB/IgHKLJ.fa" %SCRIPT_FOLDER

DH_DB   = "%s/germDB/IgHD.fa"  %SCRIPT_FOLDER

CH_DB	= "%s/germDB/IgHC_CH1.fa" %SCRIPT_FOLDER

dict_vgerm_db = {
	'H'   : VH_DB,
	'K'   : VK_DB,
	'L'   : VL_DB,
	'KL'  : VKL_DB,
	'HKL' : VHKL_DB
}

dict_jgerm_db = {
	'H'   : JH_DB,
	'K'   : JK_DB,
	'L'   : JL_DB,
	'KL'  : JKL_DB,
	'HKL' : JHKL_DB
}


ALL_FOLDERS = ["work/annotate", "work/lineage", "work/internal", "work/annotate/vgene", "work/annotate/jgene", "work/lineage/last_round", "output/sequences", "output/sequences/amino_acid", "output/sequences/nucleotide", "output/tables", "output/plots", "output/logs" ]


CMD_BLAST           = "%s -db %s -query %s -out %s -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand' -gapopen 5 -gapextend 2 -penalty -1 -reward 1 -evalue 1e-3 -max_target_seqs 10 -word_size %d"
V_BLAST_WORD_SIZE   = 7
J_BLAST_WORD_SIZE   = 5
BLAST_OTHER_OPTIONS = "-evalue 1e-3"

PBS_STRING = "\
#!/bin/bash\n\
#$ -N %s		# job name\n\
#$ -l mem=%s,time=%s	# resource requests\n\
#$ -cwd				# use current directory as job status output\n\
#$ -o /dev/null			# use sane outputs for array jobs\n\
#$ -e /dev/null\n\
\n\
%s\n"


CMD_BLASTCLUST	= "/ifs/home/c2b2/bh_lab/shares/blast/current/ia32-linux/bin/blastclust -p F -L .9 -S 95 -i %s -o %s"	#pF: nucleotide; L.9: 90%[coverage]  S: Identities 


PARSED_BLAST_HEADER = ["qid", "sid", "identity", "align_len", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "score", "strand","other_sids"]
PARSED_BLAST_HEADER_VERBOSE = ["query_id", "sbjct_id", "strand", "evalue", "score", "identities", "gaps", "aln_len", 
								"query_start", "query_end", "query_len", "sbjct_start", "sbjct_end", "aln_query", "aln_sbjct"]

