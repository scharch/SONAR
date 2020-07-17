#!/usr/bin/env python3

"""
parse_blast.py

This script is called internally from 1.3-finalize_assignments.py to allow the parsing 
    to be parallelized on the cluster.

Usage:	parse_blast.py --jmotif MOTIF --nterm OPT --chunk NUM [options]

Options:
    --jmotif TGGGG
    --nterm OPT
    --chunk NUM
    --noFallBack [default: False]

Split out from original 1.3-finalize_assignments.py by Chaim A Schramm on 2019-04-01.
Added `sequence_alignment` field for noJ reads as v gene region found by BLAST by CAS 2019-05-08.
Added `locus`, `rev-comp`, and `productive` fields for noJ reads by CA Schramm 2019-05-23.
Added locus consistency checks by CAS 2020-01-02.
Added `complete_vdj` flag by CAS 2020-07-16.

Copyright (c) 2019-2020 Vaccine Research Center, National Institutes of Health, USA.
All rights reserved.

"""

import sys, os
from docopt import docopt
from collections import Counter
import airr

try:
	from SONAR.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/annotate")
	sys.path.append(find_SONAR[0])
	from SONAR.annotate import *


def find_cdr3_borders(v_id,vgene,vlength,vstart,vend,jgene,jstart,j_start_on_read,jgaps,read_sequence):

	'''
	v_id = name of assigned V gene (eg IGHV1-2*02)
	vgene = germline sequence of assigned V gene
	vlength = length of QUERY sequence taken up by match
		(might be different from blast-reported length 
		and/or vend-vstart+1 because of in-dels)
	vstart = position on germline V gene where match begins (hopefully = 1)
	vend = position on germline V gene where match ends
	jgene = germline sequence of assigned J gene
	jstart = position on germline J gene where match begins
	j_start_on_read = position on query (v-cut version, not full 454 read) 
		where match with germline J begins
	jgaps = blast-reported number of gaps in J assignment
	read_sequence = V(D)J-trimmed sequence of the 454 read
	'''

	vMatches = []
	cys_pat = "TG[T|C|N]" #N is for a couple of shorter V's, like VH4-31
	if re.match("IGLV2-(11|23)", v_id):
		cys_pat = "TGCTGC" #special case
	if re.match("IGHV1-C",v_id):
		cys_pat = "TATGC"
	for m in re.finditer(cys_pat,vgene,flags=re.I):
		vMatches.append(m)

	#last one **IN FRAME** is the cysteine we want! (matters for light chains)
	cdr3_start=-1
	vMatches.reverse()
	for cys in vMatches:
		if cys.start() % 3 == 0:
			cdr3_start = vlength - (vend - cys.start())
			break
		
	# If BLAST has truncated the V gene alignment prior to reaching the conserved cysteine, but still found the J gene,
	#	 that likely indicates a large in-del, which must be accounted for, or the start position of CDR 3 will be wrong.
	# The only easy/automatic tool we have is to look for the conserved CxR/K motif; if it's mutated (or we are looking
	#	 at a light chain), there's nothing to do. In this case, marking CDR3 as not found is probably preferable to 
	#	 keeping the uncorrected sequence, though it should be small effect either way.
	if cdr3_start > vlength:
		has_cxrk, c_start, c_end = has_pat(read_sequence[vlength:], pat_nuc_cxrk)
		if has_cxrk:
			cdr3_start = vlength + c_start
		else:
			cdr3_start = -1

	jMatch = re.search(arguments['--jmotif'],jgene,flags=re.I)
	WF_motif = -1 #pass back to main program to check for out-of-frame junctions

	try:
		cdr3_end = vlength + j_start_on_read + (jMatch.start() - jstart) +3
		WF_motif = jMatch.start()
	except:
		cdr3_end = -1 #if we didn't find the motif, we'll count it as a bad cdr3 without crashing

	if jgaps > 0:
		#check for jMotif on read to correct for gaps and get the last one
		#but only check from cdr3 start on and don't let the end move more
		#		than a codon, because there can be similar motifs in CDR3
		wgxg = []
		for m in re.finditer(arguments['--jmotif'], read_sequence[cdr3_start:]):
			wgxg.append(m)
		if len(wgxg) > 0:
			if abs(wgxg[-1].start() + 3 - cdr3_end) <= 3:
				cdr3_end = wgxg[-1].start() + 3

	return cdr3_start, cdr3_end, WF_motif

def main():

	print( "Processing chunk %s..." % arguments['--chunk'])

	#get raw seq stats from temp table
	raw = csv.reader(open("%s/lookup_%s.txt" % (prj_tree.internal, arguments['--chunk']),'r'), delimiter=sep)

	raw_count, total, found, noV, noJ, f_ind = 0, 0, 0, 0, 0, 1
	counts = Counter()

	writer = csv.writer(open("%s/jtophit_%s.txt" %(prj_tree.jgene, arguments['--chunk']), "w"), delimiter = sep, dialect='unix', quoting=csv.QUOTE_NONE)
	writer.writerow(PARSED_BLAST_HEADER)
	dict_jcounts = dict()
	dict_ccounts = dict()
	dict_dcounts = dict()
		
	c = False
	if os.path.isfile("%s/%s_C_%s.txt" % (prj_tree.jgene, prj_name, arguments['--chunk'])):
		c = True
		cWriter = csv.writer(open("%s/ctophit_%s.txt" %(prj_tree.jgene, arguments['--chunk']), "w"), delimiter = sep, dialect='unix', quoting=csv.QUOTE_NONE)
		cWriter.writerow(PARSED_BLAST_HEADER)

	d = False
	if os.path.isfile("%s/%s_D_%s.txt" % (prj_tree.jgene, prj_name, arguments['--chunk'])):
		d = True
		dWriter = csv.writer(open("%s/dtophit_%s.txt" %(prj_tree.jgene, arguments['--chunk']), "w"), delimiter = sep, dialect='unix', quoting=csv.QUOTE_NONE)
		dWriter.writerow(PARSED_BLAST_HEADER)


	seq_stats = airr.create_rearrangement( "%s/rearrangements_%s.tsv"%(prj_tree.internal, arguments['--chunk']), fields=['complete_vdj','vj_in_frame','stop_codon','locus','c_call','junction_length','source_file','source_id','duplicate_count','length_raw','length_trimmed','indels','status','blast_identity','consensus_count','cell_id'])

	dict_vgerm_aln, dict_other_vgerms, dict_vcounts = get_top_hits("%s/%s_%s.txt"%(prj_tree.vgene, prj_name, arguments['--chunk']) )
	dict_jgerm_aln, dict_other_jgerms, dict_jcounts = get_top_hits("%s/%s_%s.txt"%(prj_tree.jgene, prj_name, arguments['--chunk']), topHitWriter=writer, dict_germ_count=dict_jcounts, strand="plus" )

	if c:
		minCStartPos = dict( [ (x, dict_jgerm_aln[x].qend) for x in dict_jgerm_aln.keys() ] )
		dict_cgerm_aln, dict_other_cgerms, dict_ccounts = get_top_hits("%s/%s_C_%s.txt"%(prj_tree.jgene, prj_name, arguments['--chunk']), topHitWriter=cWriter, dict_germ_count=dict_ccounts, minQStart=minCStartPos, strand="plus" )

	if d:
		maxDEndPos = dict( [ (x, dict_jgerm_aln[x].qstart) for x in dict_jgerm_aln.keys() ] )
		dict_dgerm_aln, dict_other_dgerms, dict_dcounts = get_top_hits("%s/%s_D_%s.txt"%(prj_tree.jgene, prj_name, arguments['--chunk']), topHitWriter=dWriter, dict_germ_count=dict_dcounts, maxQEnd=maxDEndPos, strand="plus" )

	for entry in SeqIO.parse( "%s/%s_%s.fasta" % (prj_tree.vgene, prj_name, arguments['--chunk']), "fasta"):
		total += 1

		raw_stats = next(raw)
		raw_count += 1
			
		while not entry.id == raw_stats[0]:
			#we found a read that did not meet the length cut-off
			raw_stats = next(raw)
			raw_count += 1

				
		rearrangement = dict()
		rearrangement['sequence_id'] = raw_stats[0]
		rearrangement['source_file'] = raw_stats[1]
		rearrangement['source_id']   = raw_stats[2]
		rearrangement['length_raw']  = raw_stats[3]
		rearrangement['sequence']    = str(entry.seq)
		rearrangement['complete_vdj']= False

		if not raw_stats[4] == "NA":
			rearrangement['duplicate_count'] = raw_stats[4]
		if not raw_stats[5] == "NA":
			rearrangement['consensus_count'] = raw_stats[5]
		if not raw_stats[6] == "NA":
			rearrangement['cell_id'] = raw_stats[6]
				
		if not entry.id in dict_vgerm_aln:
			noV+=1
			rearrangement['status'] = 'noV'
			seq_stats.write(rearrangement)
		elif not entry.id in dict_jgerm_aln:
			noJ+=1
			myV = dict_vgerm_aln[entry.id]
			entry.seq = entry.seq[ myV.qstart - 1 : myV.qend ]
			if (myV.strand == 'minus'):
				entry.seq = entry.seq.reverse_complement()
				rearrangement['rev_comp']       = "T"
			else:
				rearrangement['rev_comp']       = "F"
			myVgenes = ",".join( [myV.sid] + dict_other_vgerms.get(entry.id,[]) )
			
			vlocus = ""
			if re.search( "(HV|VH|heavy)", myV.sid, re.I ):
				vlocus = "IGH"
			elif re.search( "(LV|VL|lambda)", myV.sid, re.I ):
				vlocus = "IGL"
			elif re.search( "(KV|VK|kappa)", myV.sid, re.I ):
				vlocus = "IGK"

			rearrangement['v_call'] = myVgenes
			rearrangement['locus']  = vlocus
			rearrangement['productive'] = "F"
			rearrangement['status'] = 'noJ'
			rearrangement['sequence_alignment'] = str(entry.seq)
			seq_stats.write(rearrangement)

		else:
				
			found += 1
			myV = dict_vgerm_aln[entry.id]
			myJ = dict_jgerm_aln[entry.id]
			added5 = 0
			productive = "T"
			indel = "F"
			stop = "F"
			cdr3 = True
			
			vlocus = ""
			if re.search( "(HV|VH|heavy)", myV.sid, re.I ):
				vlocus = "IGH"
			elif re.search( "(LV|VL|lambda)", myV.sid, re.I ):
				vlocus = "IGL"
			elif re.search( "(KV|VK|kappa)", myV.sid, re.I ):
				vlocus = "IGK"

			#get actual V(D)J sequence
			v_len = myV.qend - (myV.qstart-1) #need to use qstart and qend instead of alignment to account for gaps

			#try to recover 3' of J
			if myJ.send < len(dict_j[myJ.sid].seq) and \
				 ( (myV.strand == "plus" and myV.qstart + v_len + myJ.qend + (len(dict_j[myJ.sid].seq)-myJ.send) <= len(entry.seq)) or \
					(myV.strand == "minus" and myV.qend - (v_len + myJ.qend + (len(dict_j[myJ.sid].seq)-myJ.send)) >= 0) ):
					vdj_len = v_len + myJ.qend + (len(dict_j[myJ.sid].seq) - myJ.send)
			else:
				vdj_len = v_len + myJ.qend

			#check for complete VDJ
			if min(myV.sstart, myV.send) == 1 and max(myJ.sstart, myJ.send) >= len(dict_j[myJ.sid].seq)-1: #-1 because the last nucleotide is part of the constant region
				rearrangement['complete_vdj'] = True

			const_seq = ""
			if (myV.strand == 'plus'):
				const_seq = str( entry.seq[myV.qstart+vdj_len-1 : ] )
				if myV.sstart > 1:
					if arguments['--nterm'] == "extend":
						if myV.qstart >= myV.sstart:
							entry.seq = entry.seq[ myV.qstart - myV.sstart : myV.qstart + vdj_len - 1 ]
							added5 = myV.sstart - 1
						else:
							entry.seq = entry.seq[  : myV.qstart + vdj_len - 1 ]
							added5 = myV.qstart - 1
					elif arguments['--nterm'] == "germline":
						entry.seq = dict_v[myV.sid].seq[ 0 : myV.sstart-1 ] + entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]
						added5 = myV.sstart - 1
					else:
								entry.seq = entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]

				else: #blast found full V gene
					entry.seq = entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]

			else: #minus strand
				const_seq = str( entry.seq[ : myV.qend-vdj_len].reverse_complement() )
				if myV.send > 1:
					if arguments['--nterm'] == "extend":
						if len(entry.seq)-myV.qend >= myV.send-1:
							entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend+myV.send-1 ].reverse_complement()
							added5 = myV.send - 1
						else:
							added5 = len(entry.seq) - myV.qend
							entry.seq = entry.seq[ myV.qend - vdj_len :  ].reverse_complement()
					elif arguments['--nterm'] == "germline":
						entry.seq = dict_v[myV.sid].seq[ 0 : myV.send-1 ] + entry.seq[ myV.qend - vdj_len : myV.qend ].reverse_complement()
						added5 = myV.send - 1
					else:
						entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend ].reverse_complement()

				else: #blast found full V gene
					entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend ].reverse_complement()

			#get CDR3 boundaries
			cdr3_start,cdr3_end,WF_motif = find_cdr3_borders(myV.sid,str(dict_v[myV.sid].seq), v_len, min(myV.sstart, myV.send), max(myV.sstart, myV.send), str(dict_j[myJ.sid].seq), myJ.sstart, myJ.qstart, myJ.gaps, str(entry.seq[ added5 : ])) #min and max statments take care of switching possible minus strand hit
			cdr3_seq = entry.seq[ added5+cdr3_start : added5+cdr3_end ]

			#push the sequence into frame for translation, if need be
			v_frame = ( min([myV.sstart, myV.send]) - added5 ) % 3
			five_prime_add = (v_frame-1) % 3
			entry.seq = 'N' * five_prime_add + entry.seq 

			#prevent BioPython errors by trimming to last full codon
			#if (len(entry.seq) % 3) > 0:
			# entry.seq = entry.seq [ : -1 * (len(entry.seq) % 3) ]

			#check for stop codons
			if '*' in entry.seq.translate():
				stop = "T"

			#check for in-frame junction
			if len(cdr3_seq) % 3 != 0:
				productive = "F"
			else: #even if recombination looks ok, might be (sequencing) indels in V and/or J
				j_frame = 3 - ( ( WF_motif - myJ.sstart ) % 3 ) #j genes start in different frames, so calculate based on position of conserved W/F found by the cdr3 subroutine above
				frame_shift = (v_len + myJ.qstart + added5 - 1) % 3

				if (v_frame + frame_shift) % 3 != j_frame % 3:
					indel = "T" 
				else:
					#use blast gaps to detect frame shift in-dels
					#most of these have stop codons or other sequence problems, but we'll catch a few extra this way
					if (abs(myV.send-myV.sstart)-(myV.qend-myV.qstart)) % 3 != 0 or ((myJ.send-myJ.sstart)-(myJ.qend-myJ.qstart)) % 3 != 0:
						indel = "T"

			#make sure cdr3 boundaries make sense
			if (cdr3_end<=cdr3_start or cdr3_end>vdj_len or cdr3_start<0):
				cdr3 = False

			status = "good"
			if not cdr3:
				status = "noCDR3"
			elif productive == "F":
				status = "nonproductive"
			elif indel == "T":
				status = "indel"
			elif stop == "T":
				status = "stop"
			elif arguments['--nterm'] == "discard" and min(myV.sstart,myV.send) > 1:
				status = "missingNterm"

			#add germline assignments to fasta description and write to disk
			myVgenes = ",".join( [myV.sid] + dict_other_vgerms.get(entry.id,[]) )
			myJgenes = ",".join( [myJ.sid] + dict_other_jgerms.get(entry.id,[]) )
				
			myDgenes = ""
			if d:
				if entry.id in dict_dgerm_aln:
					if not vlocus in ["IGK", "IGL"]:
						#supress spurious D gene hits if it's a light chain
						myDgenes = ",".join( [dict_dgerm_aln[entry.id].sid] + dict_other_dgerms.get(entry.id,[]) )

			myCgenes = ""
			if c and entry.id in dict_cgerm_aln:
				myCgenes = ",".join( [dict_cgerm_aln[entry.id].sid] + dict_other_cgerms.get(entry.id,[]) )
			elif not arguments['--noFallBack']:
				if re.match("C[CT]", const_seq):
					myCgenes = "IGHG" #could also be IgE, but I'm assuming that's rare
				elif re.match("GGA", const_seq):
					myCgenes = "IGHM"
				elif re.match("CAT", const_seq):
					myCgenes = "IGHA"
				elif re.match("CAC", const_seq):
					myCgenes = "IGHD"
				elif re.match("CGA", const_seq):
					myCgenes = "IGKC"
				elif re.match("GGT", const_seq):
					myCgenes = "IGLC"

			jlocus = ""
			if re.search( "(HJ|JH|heavy)", myJ.sid, re.I ):
				jlocus = "IGH"
			elif re.search( "(LJ|Jl|lambda)", myJ.sid, re.I ):
				jlocus = "IGL"
			elif re.search( "(KJ|JK|kappa)", myJ.sid, re.I ):
				jlocus = "IGK"

			if not vlocus == jlocus:
				#this really shouldn't happen unless one or both gene assignments are
				#    based on very short partial hits. Unfortuantely, the lengths/e-values
				#    are on different scales, so I don't currently have a good heuristic to
				#    pick between the two. Just flag it and give up, at least for now.
				status = "chimera"

			if not myCgenes == "" and not vlocus in myCgenes: #will fail for custom libraries where C gene names don't start with locus
				myCgenes = "" #assume constant is incorrect since usually based on only a few bases

			#do AIRR output
			if myV.strand == "plus":
				rearrangement['rev_comp']       = "F"
			else:
				rearrangement['rev_comp']       = "T"
			if status == "good":
				rearrangement['productive']     = "T"
			else:
				rearrangement['productive']     = "F"
			rearrangement['vj_in_frame']        = productive
			rearrangement['stop_codon']         = stop
			rearrangement['locus']              = vlocus
			rearrangement['v_call']             = myVgenes
			rearrangement['j_call']             = myJgenes
			rearrangement['d_call']             = myDgenes
			rearrangement['c_call']             = myCgenes
			rearrangement['sequence_alignment'] = str(entry.seq)
			rearrangement['junction']           = cdr3_seq
			rearrangement['junction_aa']        = cdr3_seq.translate()
			rearrangement['junction_length']    = len(cdr3_seq)
			rearrangement['length_trimmed']     = len(entry.seq)
			rearrangement['indels']             = indel
			rearrangement['status']             = status
			rearrangement['blast_identity']     = "%.3f" % (myV.identity/100.0)
				
			seq_stats.write(rearrangement)

			counts[status] += 1

	print( "chunk %s: %d done, found %d; %d good..." %(arguments['--chunk'], total, found, counts['good']) )

	seq_stats.close()


if __name__ == '__main__':
	
	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	#load saved locus and library information
	handle  = open( "%s/gene_locus.txt" % prj_tree.internal, "r")
	species = handle.readline().strip()
	locus   = handle.readline().strip()
	vlib    = handle.readline().strip()
	jlib    = handle.readline().strip()

	dict_v = load_fastas(vlib)
	dict_j = load_fastas(jlib)

	main()
