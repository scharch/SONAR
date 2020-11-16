#!/usr/bin/env python3


"""
Core functions for SONAR

Created by Zhenhai Zhang on 2011-04-05 as mytools.py
Edited and commented for publication by Chaim A Schramm on 2015-02-10.
Added column comparisons and accept a RearrangementReader object for
     filterAirrTsv by CA Schramm on 2020-06-11.
Changed filterAirrTsv to use eval by CA Schramm on 2020-07-02.

Copyright (c) 2011-2020 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.
"""

from Bio import SeqIO
from Bio import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

import airr

from math import log
from io import StringIO

import glob
import re
import subprocess
import atexit


from numpy import mean, array, zeros, ones, nan, std, isnan

from SONAR.commonVars import *


########## COMMAND LOGGING ############
global printLog
printLog = False #this tells us whether or not to print log info on exit (skip if program was called with -h)


#overload default error handling so we can log whether it was a successful exit or not
# code taken from: https://stackoverflow.com/a/9741784
class ExitHooks(object):
    def __init__(self):
        self.exit_code = None
        self.exception = None

    def hook(self):
        self._orig_exit = sys.exit
        sys.exit = self.exit
        self._orig_except = sys.excepthook
        sys.excepthook = self.exc_handler

    def exit(self, code=0):
        self.exit_code = str(code)
        self._orig_exit(code)

    def exc_handler(self, exc_type, exc, tb):
        self.exception = exc_type.__name__ + ": " + str(exc)
        self._orig_except(exc_type, exc, tb)

hooks = ExitHooks()
hooks.hook()


def logCmdLine( command ):

    global printLog, logFile

    logFile = "SONAR_command_history.log"
    if os.path.isdir( "%s/output/logs" % os.getcwd() ):
        logFile = "%s/output/logs/command_history.log"%os.getcwd()

    for idx,arg in enumerate(command):
        if re.search("(\s|\*)", arg):
            command[idx] = "'"+arg+"'"

    p = subprocess.Popen(['git', '--git-dir', os.path.dirname(command[0])+"/../.git",
    					  '--work-tree', os.path.dirname(command[0])+"/../",
                          'describe', '--always','--dirty','--tags'],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    VERSION = p.communicate()[0].decode().strip()

    logStatement = "\n%s -- SONAR %s run with command:\n\t%s\n" % (time.strftime("%c"), VERSION, " ".join(command))
    print(logStatement, file=sys.stderr)

    try:
        with open(logFile, "a") as handle:
            handle.write( logStatement )

        printLog = True

    except:
        print("Directory appears to be read-only; command line and output will not be saved", file=sys.stderr)


def logExit():

    global printLog, logFile
    if printLog:
        with open(logFile, "a") as handle:
            if hooks.exit_code is not None:
                formatted = re.sub( "\n", "\n\t", hooks.exit_code.strip(" \t\r\n") ) #remove white space and new lines on both ends; indent if multiple lines
                handle.write( "%s -- Program exited with error:\n\t%s\n" % (time.strftime("%c"),formatted) )
            elif hooks.exception is not None:
                formatted = re.sub( "\n", "\n\t", hooks.exception.strip(" \t\r\n") ) #remove white space and new lines on both ends; indent if multiple lines
                handle.write( "%s -- Exception:\n\t%s\n" % (time.strftime("%c"),formatted) )
            else:
                handle.write( "%s -- Program finished successfully\n" % time.strftime("%c") )

atexit.register(logExit)




#
# -- BEGIN -- class definition
#

class MyAlignment:
	def __init__(self, row):
		self.qid	= row[0].strip()		# query id
		self.sid	= row[1].strip()		# subject id
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



class ProjectFolders:
	"""folder structure of a project """

	def __init__(self, proj_home, create=True):

		self.home       = proj_home
		self.work       = "%s/work"        %  proj_home
		self.out        = "%s/output"      %  proj_home

		#working folders
		self.preprocess = "%s/preprocess"  %  self.work
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

		if create:
		    os.makedirs( self.preprocess, exist_ok=True )
		    os.makedirs( self.internal,   exist_ok=True )
		    os.makedirs( self.vgene,      exist_ok=True )
		    os.makedirs( self.jgene,      exist_ok=True )
		    os.makedirs( self.last,       exist_ok=True )
		    os.makedirs( self.beast,      exist_ok=True )
		    os.makedirs( self.aa,         exist_ok=True )
		    os.makedirs( self.nt,         exist_ok=True )
		    os.makedirs( self.tables,     exist_ok=True )
		    os.makedirs( self.plots,      exist_ok=True )
		    os.makedirs( self.logs,       exist_ok=True )
		    os.makedirs( self.rates,      exist_ok=True )
#
# -- END -- class definition
#


#
# -- BEGIN --  folder and file methods
#

def fullpath2last_folder(s):
	"""get immediate parent folder"""

	return s[s.rindex("/") + 1 :]


#
# -- END -- folder and file methods
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


def load_seqs_in_dict(f, ids):
	"""
	load all sequences in file f is their id is in ids (list or set or dictionary)
	"""
	result = dict()
	for entry in SeqIO.parse(open(f, "r"), "fasta"):
		if entry.id in ids:
			result[entry.id] = entry

	return result


def load_fastas_in_list(f, l):

	print( "loading reads from %s as in given list..." %f )
	reader, result, good = SeqIO.parse(open(f, "r"), "fasta"), dict(), 0

	for entry in reader:
		if entry.id in l:

			#changed to match load from set CAS 20121004
			result[entry.id] = entry
			good += 1
			if good == len(l): break

			#myseq = MySeq(entry.id, entry.seq)
			#myseq.desc = entry.description
			#result[entry.id] = myseq

	print( "%d loaded...." %len(result) )
	return result


def load_fastas_with_Vgene(f, v):
	print( "loading reads from %s assigned to %s..." %(f,v) )
	reader, dict_reads = SeqIO.parse(open(f, "r"), "fasta"), dict()
	for entry in reader:
		if re.search(v, entry.description):
			dict_reads[entry.id] = entry

	print( "%d loaded..." %len(dict_reads) )
	return dict_reads


def load_fastas(f):
	"""return gene ID and sequences in a dictionary"""
	print( "loading sequence info from %s..." %f )

	reader, result = SeqIO.parse(open(f, "r"), "fasta"), dict()

	for entry in reader:
		#myseq = MySeq(entry.id, entry.seq)
		#myseq.desc = entry.description
		result[entry.id] = entry

	return result


def generate_read_fasta(f):
	"""read fasta file and yield one reads per time """

	filetype = "fasta"
	if re.search("\.(fq|fastq)$", f) is not None:
		filetype = "fastq"

	reader = SeqIO.parse(open(f, "r"), filetype)
	for entry in reader:
		yield entry


def generate_read_fasta_folder(fastas):

	for fasta_file in fastas:

		filetype = "fasta"
		if re.search("\.(fq|fastq)$", fasta_file) is not None:
			filetype = "fastq"

		for entry in SeqIO.parse(open(fasta_file, "r"), filetype):

			yield entry, None, fasta_file

			#if we reimplement qual handling uncomment next section
			#if filetype == "fastq":
			#	yield entry, MyQual(entry.id, entry.letter_annotations["phred_quality"]), fasta_file



def check_fasta_qual_pair(fa_name, qu_name):
	return fa_name[ : fa_name.rindex(".")] == qu_name[ : qu_name.rindex(".")]


def translate_a_sequence(s):
	"""translate nucleotide to protein"""  # this method is redundant, but we wannt deal with "N"
	s = s.upper()
	if s.find("N") > 0:
		return None
	codons = ["".join(x) for x in zip(s[ :: 3], s[1 :: 3], s[2 :: 3])]

	return "".join([dict_codon2aa[x][0] for x in codons])

#
# -- END -- FASTA file and sequence methods
#


#
# -- BEGIN -- alignment functions
#

def quickAlign( refseq, testseq, maxiters=None, diags=None, gapopen=None ):

	#sanity check
	try:
		refseq	= re.sub( "-", "", refseq )
	except TypeError:
		#not a string, probably a SeqRecord
		try:
			refseq = str( refseq.seq )
			refseq	= re.sub( "-", "", refseq )
		except AttributeError:
			#give up
			sys.exit( "quickAlign() requires inputs to be either strings or SeqRecord objects" )

	try:
		testseq	= re.sub( "-", "", testseq )
	except TypeError:
		#not a string, probably a SeqRecord
		try:
			testseq = str( testseq.seq )
			testseq	= re.sub( "-", "", testseq )
		except AttributeError:
			#give up
			sys.exit( "quickAlign() requires inputs to be either strings or SeqRecord objects" )

	handle = StringIO()
	handle.write( ">ref\n%s\n>test\n%s\n"%(refseq,testseq) )
	data = handle.getvalue()

	muscle_cline = MuscleCommandline(cmd=muscle, quiet=True)
	if maxiters is not None: muscle_cline.maxiters = maxiters
	if diags    is not None: muscle_cline.diags    = diag
	if gapopen  is not None: muscle_cline.gapopen  = gapopen

	stdout, stderr = muscle_cline(stdin=data)

	aligned = dict()
	for p in SeqIO.parse(StringIO(stdout), "fasta"):
		aligned[ p.id ] = str(p.seq)
	return aligned


def scoreAlign( alignDict, reference="ref", query="test", countTerminalGaps=False, countInternalGaps=True, skip=0 ):

	refLen = len( re.sub("-", "", alignDict[reference]) )

	if not countTerminalGaps:
		leftGap = re.match( "-+", alignDict[reference] )
		if not leftGap:
			leftGap = re.match( "-+", alignDict[query] )
		if leftGap:
			alignDict = { s:alignDict[s][leftGap.end():] for s in [reference,query] }
		rightGap = re.search( "-+$", alignDict[reference] )
		if not rightGap:
			rightGap = re.search( "-+$", alignDict[query] )
		if rightGap:
			alignDict = { s:alignDict[s][0:rightGap.start()] for s in [reference,query] }

	position = 0
	alignLen = 0.0
	covNum   = 0.0
	matches	 = 0
	for r,t in zip(alignDict[reference], alignDict[query]):
		if not (r == "-" or t == "-"):
			covNum += 1
		if (not countInternalGaps) and (r == "-" or t == "-"):
			continue
		elif position < skip:
			position += 1
			continue
		else:
			alignLen += 1
			if r == t:
				matches += 1

		coverage = covNum  / refLen

	if alignLen == 0:
	    return 0, 0
	else:
	    return matches/alignLen, coverage

#
# -- END -- alignment functions
#


#
# -- BEGIN -- AIRR manipulation functions
#

def filterAirrTsv(rearrangementsFile, ruleList, useOR=False):
	good = 0

	try:
		#see if it's a file name
		reader = airr.read_rearrangement( rearrangementsFile )
	except TypeError:
		#assume it's an already open RearrangementReader object instead
		reader = rearrangementsFile

	for r in reader:
		keep = False
		if useOR:
			keep = any( [eval(rule, {'re':re}, {'r':r}) for rule in ruleList] )
		else:
			keep = all( [eval(rule, {'re':re}, {'r':r}) for rule in ruleList] )

		if keep:
			good += 1
			if good % 10000 == 0:
				sys.stderr.write("Found %d matching rearrangements so far...\n" % good)
			yield r


def airrToFasta( rearrangements, field='sequence_alignment', aa=False):
	for r in rearrangements:
		if r[field] == "":
			continue

		tempSeq = SeqRecord( id=r['sequence_id'], seq=Seq.Seq(re.sub("[-.+]","",r[field])) )
		if aa:
			tempSeq.seq = tempSeq.seq.translate()

		def_line = ""
		if not r['v_call'] == '':          def_line += " v_call=%s"          % r['v_call']
		if not r['d_call'] == '':          def_line += " d_call=%s"          % r['d_call']
		if not r['j_call'] == '':          def_line += " j_call=%s"          % r['j_call']
		if not r['junction'] == '':        def_line += " junction=%s"        % r['junction']
		if not r['junction_aa'] == '':     def_line += " junction_aa=%s"     % r['junction_aa']
		if not r['junction_length'] == '': def_line += " junction_length=%s" % r['junction_length']
		if 'locus' in r and not r['locus']  == '':                    def_line += " locus=%s"           % r['locus']
		if 'c_call' in r and not r['c_call'] == '':                   def_line += " c_call=%s"          % r['c_call']
		if 'status' in r and not r['status'] == '':                   def_line += " status=%s"          % r['status']
		if 'v_identity' in r and not r['v_identity'] == '':           def_line += " v_identity=%s"      % r['v_identity']
		if 'duplicate_count' in r and not r['duplicate_count'] == '': def_line += " duplicate_count=%s" % r['duplicate_count']
		if 'consensus_count' in r and not r['consensus_count'] == '': def_line += " consensus_count=%s" % r['consensus_count']
		if 'cell_id' in r and not r['cell_id'] == '':                 def_line += " cell_id=%s"         % r['cell_id']
		if 'cell_status' in r and not r['cell_status'] == '':         def_line += " cell_status=%s"     % r['cell_status']
		if 'centroid' in r and not r['centroid'] == '':               def_line += " centroid=%s"        % r['centroid']
		if 'cluster_count' in r and not r['cluster_count'] == '':     def_line += " cluster_count=%s"   % r['cluster_count']
		if 'clone_id' in r and not r['clone_id'] == '':               def_line += " clone_id=%s"        % r['clone_id']
		if 'clone_count' in r and not r['clone_count'] == '':         def_line += " clone_count=%s"     % r['clone_count']
		tempSeq.description = def_line

		yield tempSeq

#
# -- END -- AIRR manipulation functions
#
