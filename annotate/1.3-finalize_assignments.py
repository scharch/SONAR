#!/usr/bin/env python3

"""
1.3-finalize_assignments.py

This script parses the BLAST output from 1.1-blast-V_assignment.py and 
      1.2-blast_J.py. Sequences with successful assignments are
      output into fasta files and a master table is created summarizing the
      properties of all input sequences.

Usage:  1.3-finalize_assignments.py [ --jmotif "TT[C|T][G|A]G" --nterm truncate --noclean ]

Options:
    --jmotif TGGGG    Conserved nucleotide sequence indicating the start of FWR4 on
                         the J gene. Defaults to either TGGGG for heavy chains or
                         TT[C|T][G|A]G for light chains; set manually for custom 
                         light chain libraries or for species with a different motif.
    --nterm OPT       What to do if blast hit does not extend to the N terminus of the
                         germline V gene. Options are:
                            truncate - trim to blast hit;
                            extend   - change trim boundary to correspond to expected N
                                      N terminus of germline or beginning of read if
                                      shorter (NOT recommended, typically results in
                                      bad sequences)
                            germline - replace missing region with the germline V
                                      sequence (useful for FWR1 primers)
                            discard  - only mark as 'good' sequences with full germline
                                      region
                         [default: truncate]
    --noclean         Disable automatic deletion of working files from 1.1 and 1.2 [default: False]

Created by Chaim A Schramm on 2013-07-05
Edited and commented for publication by Chaim A Schramm on 2015-02-25.
Edited to add custom J motif option for other species by CAS 2016-05-16.
Edited to add options to handle missing N terminal by CAS 2017-05-08.
Edited to make motif matching case-insensitive and to disable auto-clean-up
                               by CAS 2017-07-31.
Edited to distinguish between nonproductive junctions and apparent frameshifts
                               within V by CAS 2018-07-23.
Edited to use Py3 and DocOpt by CAS 2018-08-22.
Edited to use AIRR datarep convention by CAS 2018-10-12.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""

import sys, os
from docopt import docopt
import airr

try:
	from sonar.annotate import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/annotate")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *

global jMotif

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
	#   that likely indicates a large in-del, which must be accounted for, or the start position of CDR 3 will be wrong.
	# The only easy/automatic tool we have is to look for the conserved CxR/K motif; if it's mutated (or we are looking
	#   at a light chain), there's nothing to do. In this case, marking CDR3 as not found is probably preferable to 
	#   keeping the uncorrected sequence, though it should be small effect either way.
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
		#    than a codon, because there can be similar motifs in CDR3
		wgxg = []
		for m in re.finditer(arguments['--jmotif'], read_sequence[cdr3_start:]):
			wgxg.append(m)
		if len(wgxg) > 0:
			if abs(wgxg[-1].start() + 3 - cdr3_end) <= 3:
				cdr3_end = wgxg[-1].start() + 3

	return cdr3_start, cdr3_end, WF_motif


def main():

	if not glob.glob("%s/%s_*.fasta" % (prj_tree.jgene, prj_name)):
		sys.exit("No jBlast output found!\n")
	
	print( "curating junction and 3' end..." )


	allV_aa	     = open ("%s/%s_allV.fa"	 % (prj_tree.aa, prj_name), "w" )
	allV_nt	     = open( "%s/%s_allV.fa"	 % (prj_tree.nt, prj_name), "w" )

	allJ_aa	     = open( "%s/%s_allJ.fa"	 % (prj_tree.aa, prj_name), "w" )
	allJ_nt	     = open( "%s/%s_allJ.fa"	 % (prj_tree.nt, prj_name), "w" )

	vj_aa	     = open( "%s/%s_goodVJ.fa"	 % (prj_tree.aa, prj_name), "w" )
	vj_nt	     = open( "%s/%s_goodVJ.fa"	 % (prj_tree.nt, prj_name), "w" )

	good_cdr3_aa = open( "%s/%s_goodCDR3.fa" % (prj_tree.aa, prj_name), "w" )
	good_cdr3_nt = open( "%s/%s_goodCDR3.fa" % (prj_tree.nt, prj_name), "w" )

	all_cdr3_nt  = open( "%s/%s_allCDR3.fa"	 % (prj_tree.nt, prj_name), "w" )


	#get raw seq stats from temp table
	raw = csv.reader(open("%s/id_lookup.txt" % prj_tree.internal,'rU'), delimiter=sep)


	raw_count, total, found, noV, noJ, f_ind  = 0, 0, 0, 0, 0, 1
	counts = {'good':0,'nonproductive':0,'indel':0,'noCDR3':0,'stop':0}
	if arguments['--nterm'] == "discard":
		counts["missingNterm"]=0

	writer = csv.writer(open("%s/%s_jgerm_tophit.txt" %(prj_tree.tables, prj_name), "w"), delimiter = sep)
	writer.writerow(PARSED_BLAST_HEADER)
	dict_jcounts = dict()
	dict_ccounts = dict()
	dict_dcounts = dict()
		
	c = False
	if os.path.isfile("%s/%s_C_001.txt" % (prj_tree.jgene, prj_name)):
		c = True
		cWriter = csv.writer(open("%s/%s_cgerm_tophit.txt" %(prj_tree.tables, prj_name), "w"), delimiter = sep)
		cWriter.writerow(PARSED_BLAST_HEADER)

	d = False
	if os.path.isfile("%s/%s_D_001.txt" % (prj_tree.jgene, prj_name)):
		d = True
		dWriter = csv.writer(open("%s/%s_dgerm_tophit.txt" %(prj_tree.tables, prj_name), "w"), delimiter = sep)
		dWriter.writerow(PARSED_BLAST_HEADER)


	seq_stats = airr.create_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name), fields=['vj_in_frame','stop_codon','locus','c_call','junction_length','source_file','source_id','duplicate_count','length_raw','length_trimmed','indels','status','blast_identity'])

	
	while os.path.isfile("%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind)):

		dict_vgerm_aln, dict_other_vgerms, dict_vcounts	 =  get_top_hits("%s/%s_%03d.txt"%(prj_tree.vgene, prj_name, f_ind) )
		dict_jgerm_aln, dict_other_jgerms, dict_jcounts	 =  get_top_hits("%s/%s_%03d.txt"%(prj_tree.jgene, prj_name, f_ind), topHitWriter=writer, dict_germ_count=dict_jcounts )

		if c:
			minCStartPos = dict( [ (x, dict_jgerm_aln[x].qend) for x in dict_jgerm_aln.keys() ] )
			dict_cgerm_aln, dict_other_cgerms, dict_ccounts	 =  get_top_hits("%s/%s_C_%03d.txt"%(prj_tree.jgene, prj_name, f_ind), topHitWriter=cWriter, dict_germ_count=dict_ccounts, minQStart=minCStartPos )

		if d:
			maxDEndPos = dict( [ (x, dict_jgerm_aln[x].qstart) for x in dict_jgerm_aln.keys() ] )
			dict_dgerm_aln, dict_other_dgerms, dict_dcounts	 =  get_top_hits("%s/%s_D_%03d.txt"%(prj_tree.jgene, prj_name, f_ind), topHitWriter=dWriter, dict_germ_count=dict_dcounts, maxQEnd=maxDEndPos )

		for entry in SeqIO.parse( "%s/%s_%03d.fasta" % (prj_tree.vgene, prj_name, f_ind), "fasta"):
			total += 1

			raw_stats = next(raw)
			raw_count += 1
			
			while not entry.id == raw_stats[0]:
				#we found a read that did not meet the length cut-off
				raw_stats = next(raw)
				raw_count += 1

				
			rearrangement = dict()
			rearrangement['sequence_id']     = raw_stats[0]
			rearrangement['source_file']     = raw_stats[1]
			rearrangement['source_id']       = raw_stats[2]
			rearrangement['length_raw']      = raw_stats[4]
			rearrangement['sequence']        = str(entry.seq)

                        if not raw_stats[3] == "NA":
			        rearrangement['duplicate_count'] = raw_stats[3]
                                entry.description = "duplicate_count=%s" % raw_stats[3]
                        else:
                                entry.description = "" #just in case
				
			if not entry.id in dict_vgerm_aln:
				noV+=1
				rearrangement['status']	= 'noV'
				seq_stats.write(rearrangement)
			elif not entry.id in dict_jgerm_aln:
				noJ+=1
				myV = dict_vgerm_aln[entry.id]
				if (myV.strand == 'plus'):
					entry.seq = entry.seq[ myV.qstart - 1 :	 ]							
				else:
					entry.seq = entry.seq[	: myV.qend ].reverse_complement()
				myVgenes = ",".join( [myV.sid] + dict_other_vgerms.get(entry.id,[]) )
				entry.description += " v_call=%s status=noJ" % (myVgenes)
				allV_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq))

				#prevent BioPython errors
				if (len(entry.seq) % 3) > 0:
					entry.seq = entry.seq [ :  -1 * (len(entry.seq) % 3) ]
				allV_aa.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq.translate()))

				rearrangement['v_call'] = myVgenes
				rearrangement['status']	= 'noJ'
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
				
				#get actual V(D)J sequence
				v_len	= myV.qend - (myV.qstart-1) #need to use qstart and qend instead of alignment to account for gaps

				#try to recover 3' of J
				if myJ.send < len(dict_j[myJ.sid].seq) and \
				   ( (myV.strand == "plus" and myV.qstart + v_len + myJ.qend + (len(dict_j[myJ.sid].seq)-myJ.send) <= len(entry.seq)) or \
				     (myV.strand == "minus" and myV.qend - (v_len + myJ.qend + (len(dict_j[myJ.sid].seq)-myJ.send)) >= 0) ):
						vdj_len = v_len + myJ.qend + (len(dict_j[myJ.sid].seq) - myJ.send)
				else:
					vdj_len = v_len + myJ.qend
					
				if (myV.strand == 'plus'):
					if myV.sstart > 1:
						if arguments['--nterm'] == "extend":
							if myV.qstart >= myV.sstart:
								entry.seq = entry.seq[ myV.qstart - myV.sstart : myV.qstart + vdj_len - 1 ]
								added5 = myV.sstart - 1
							else:
								entry.seq = entry.seq[	: myV.qstart + vdj_len - 1 ]
								added5 = myV.qstart - 1
						elif arguments['--nterm'] == "germline":
							entry.seq = dict_v[myV.sid].seq[ 0 : myV.sstart-1 ] + entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]
							added5 = myV.sstart - 1
						else:
						      entry.seq = entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]

					else: #blast found full V gene
						entry.seq = entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]

				else: #minus strand
					if myV.send > 1:
						if arguments['--nterm'] == "extend":
							if len(entry.seq)-myV.qend >= myV.send-1:
								entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend+myV.send-1 ].reverse_complement()
								added5 = myV.send - 1
							else:
								added5 = len(entry.seq) - myV.qend
								entry.seq = entry.seq[	myV.qend - vdj_len :  ].reverse_complement()
						elif arguments['--nterm'] == "germline":
							entry.seq = dict_v[myV.sid].seq[ 0 : myV.send-1 ] + entry.seq[	myV.qend - vdj_len : myV.qend ].reverse_complement()
							added5 = myV.send - 1
						else:
							entry.seq = entry.seq[	myV.qend - vdj_len : myV.qend ].reverse_complement()

					else: #blast found full V gene
						entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend ].reverse_complement()

				#get CDR3 boundaries
				cdr3_start,cdr3_end,WF_motif = find_cdr3_borders(myV.sid,str(dict_v[myV.sid].seq), v_len, min(myV.sstart, myV.send), max(myV.sstart, myV.send), str(dict_j[myJ.sid].seq), myJ.sstart, myJ.qstart, myJ.gaps, str(entry.seq[ added5 :  ])) #min and max statments take care of switching possible minus strand hit
				cdr3_seq = entry.seq[ added5+cdr3_start : added5+cdr3_end ]

				#push the sequence into frame for translation, if need be
				v_frame = ( min([myV.sstart, myV.send]) - added5 ) % 3
				five_prime_add = (v_frame-1) % 3
				entry.seq = 'N' * five_prime_add + entry.seq 

				#prevent BioPython errors by trimming to last full codon
				#if (len(entry.seq) % 3) > 0:
				#	entry.seq = entry.seq [ :  -1 * (len(entry.seq) % 3) ]

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
						indel = "T"   #for gDNA we would probably want to distinguish between an out-of-frame recombination and sequencing in-dels in V or J
								#but that can be ambiguous and for cDNA we can assume that it's sll sequencing in-del anyway, even in CDR3.
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
						myDgenes = ",".join( [dict_dgerm_aln[entry.id].sid] + dict_other_dgerms.get(entry.id,[]) )

				myCgenes = ""
				if c:
					if entry.id in dict_cgerm_aln:
						myCgenes = ",".join( [dict_cgerm_aln[entry.id].sid] + dict_other_cgerms.get(entry.id,[]) )

				vlocus = ""
				if any( x in myV.sid for x in ["HV", "VH", "Vh", "vh", "heavy", "Heavy", "HEAVY"] ):
					vlocus = "IGH"
				elif any( x in myV.sid for x in ["LV", "VL", "Vl", "vl", "lambda", "Lambda", "LAMBDA"] ):
					vlocus = "IGL"
				elif any( x in myV.sid for x in ["KV", "VK", "Vk", "vk", "kappa", "Kappa", "KAPPA"] ):
					vlocus = "IGK"
					
				entry.description += " v_call=%s" % myVgenes
				if myDgenes != "":
					entry.description += " d_call=%s" % myDgenes
				entry.description += " j_call=%s" % myJgenes
				if myCgenes != "":
					entry.description += " c_call=%s" % myCgenes
				if vlocus != "":
					entry.description += " locus=%s" % vlocus
				entry.description += " status=%s blast_identity=%.3f junction_length=%d junction=%s junction_aa=%s" % ( status, myV.identity/100.0, len(cdr3_seq), cdr3_seq, cdr3_seq.translate() )

				allV_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq))
				allV_aa.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq.translate()))

				allJ_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq))
				allJ_aa.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq.translate()))

				if status == "good":

					vj_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq))
					vj_aa.write(">%s %s\n%s\n" %(entry.id, entry.description, entry.seq.translate()))

					good_cdr3_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, cdr3_seq))
					good_cdr3_aa.write(">%s %s\n%s\n" %(entry.id, entry.description, cdr3_seq.translate()))

					all_cdr3_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, cdr3_seq))

				elif cdr3:
					#CDR3 but not "good"
					all_cdr3_nt.write(">%s %s\n%s\n" %(entry.id, entry.description, cdr3_seq))


				#do AIRR output
				if myV.strand == "plus":
					rearrangement['rev_comp']   = "F"
				else:
					rearrangement['rev_comp']   = "T"
				if status == "good":
					rearrangement['productive'] = "T"
				else:
					rearrangement['productive'] = "F"
				rearrangement['vj_in_frame']	    = productive
				rearrangement['stop_codon']	    = stop
				rearrangement['locus']		    = vlocus
				rearrangement['v_call']		    = myVgenes
				rearrangement['j_call']		    = myJgenes
				rearrangement['d_call']	            = myDgenes
				rearrangement['c_call']	            = myCgenes
				rearrangement['sequence_alignment'] = str(entry.seq)
				rearrangement['junction']	    = cdr3_seq
				rearrangement['junction_aa']	    = cdr3_seq.translate()
				rearrangement['junction_length']    = len(cdr3_seq)
				rearrangement['length_trimmed']	    = len(entry.seq)
				rearrangement['indels']		    = indel
				rearrangement['status']		    = status
				rearrangement['blast_identity']	    = "%.3f" % (myV.identity/100.0)
				
				seq_stats.write(rearrangement)

				counts[status] += 1

		print( "%d done, found %d; %d good..." %(total, found, counts['good']) )
		f_ind += 1

	seq_stats.close()

	#print out some statistics
	handle = open("%s/%s_jgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
	writer	= csv.writer(handle, delimiter = sep)
	keys	= sorted(dict_jcounts.keys())
	writer.writerow(["gene", "count", "percent"])
	for key in keys:
		aline = [ key, dict_jcounts[key], "%4.2f" % (dict_jcounts[key] / float(found) * 100) ]
		writer.writerow(aline)
	handle.close()

	if len(dict_ccounts) > 0:
		handle = open("%s/%s_cgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
		writer	= csv.writer(handle, delimiter = sep)
		keys	= sorted(dict_ccounts.keys())
		writer.writerow(["gene", "count", "percent"])
		for key in keys:
			aline = [ key, dict_ccounts[key], "%4.2f" % (dict_ccounts[key] / float(found) * 100) ]
			writer.writerow(aline)
		handle.close()

	if len(dict_dcounts) > 0:
		handle = open("%s/%s_dgerm_stat.txt" %(prj_tree.tables, prj_name),'w')
		writer	= csv.writer(handle, delimiter = sep)
		keys	= sorted(dict_dcounts.keys())
		writer.writerow(["gene", "count", "percent"])
		for key in keys:
			aline = [ key, dict_dcounts[key], "%4.2f" % (dict_dcounts[key] / float(found) * 100) ]
			writer.writerow(aline)
		handle.close()

	message = "\nTotal raw reads: %d\nCorrect Length: %d\nV assigned: %d\nJ assigned: %d\nCDR3 assigned: %d\nIn-frame junction: %d\nNo indels: %d\nContinuous ORF with no stop codons: %d\n\n"  % \
							     (raw_count, total, total-noV, found, found-counts['noCDR3'], found-counts['noCDR3']-counts['nonproductive'], found-counts['noCDR3']-counts['nonproductive']-counts['indel'], counts['good'])
	print( message )
	handle = open("%s/finalize_blast.log"%prj_tree.logs, "w")
	handle.write(message)
	handle.close()

	#clean up!!
	oldFiles = glob.glob("%s/*txt"%prj_tree.vgene) + glob.glob("%s/*fasta"%prj_tree.vgene) +  glob.glob("%s/*txt"%prj_tree.jgene) + glob.glob("%s/*fasta"%prj_tree.jgene) + glob.glob("%s/id_lookup.txt"%prj_tree.internal)
	if len(oldFiles) > 0 and not arguments['--noclean']:
		[os.remove(f) for f in oldFiles]
			


if __name__ == '__main__':
	
	#log command line
	logCmdLine(sys.argv)

	arguments = docopt(__doc__)
	
	prj_tree  = ProjectFolders(os.getcwd())
	prj_name  = fullpath2last_folder(prj_tree.home)

	#load saved locus and library information
	handle = open( "%s/gene_locus.txt" % prj_tree.internal, "rU")
	locus = handle.readline().strip()
	vlib  = handle.readline().strip()
	jlib  = handle.readline().strip()

	if arguments['--jmotif'] is None:
		arguments['--jmotif'] = "TGGG"
		if "K" in locus or "L" in locus: #it's a light chain!
			if "H" in locus: #need both motifs
				arguments['--jmotif'] = "(TGGGG|TT[C|T][G|A]G)"
			else:
				arguments['--jmotif'] = "TT[C|T][G|A]G"

	if arguments['--nterm'] not in ["truncate", "extend", "germline", "discard"]:
		sys.exit("--nterm must be one of ('truncate', 'extend', 'germline', 'discard')")
	
	dict_v	  =  load_fastas(vlib)
	dict_j	  =  load_fastas(jlib)

	main()
