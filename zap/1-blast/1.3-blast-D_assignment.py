#!/usr/bin/env python
# encoding: utf-8
"""
1.1b-parse_j_blast.py

Created by Chaim A Schramm on 2013-07-05
Copyright (c) 2013 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 
1.1b-parse_j_blast.py -l <0|1|2|3>

"""
# --BEGIN-- general package
import sys
import os
import subprocess

from mytools import *
# -- END -- general package

def write_dict2file(d, total, header):
        filename = "%s/%s_jgerm_stat.txt" %(prj_tree.data, prj_name)
        writer  = csv.writer(open(filename, "w"), delimiter = sep)
        keys    = sorted(d.keys())

        writer.writerow(header)

        for key in keys:
                aline = [key, d[key], d[key] / float(total) * 100]
                writer.writerow(aline)



def get_top_hit(resultFile, germDict, writer=None):

	dict_germ_aln  =  dict()
	dict_j_seconds =  dict()

	#db_name = fullpath2last_folder(folder)

	for my_alignment, row, others in generate_blast_top_hists(resultFile):
		qid, sid 	= 	my_alignment.qid, my_alignment.sid
		aline 		= 	row + [my_alignment.strand]

		if my_alignment.strand == "-" and writer is not None:
			continue #everything was converted to positive strand when assigning v! (but we might be reprocessing V)

                aline.append(germDict[sid].seq_len)
		if len(others)>0:
			aline.append(",".join(others))
			dict_j_seconds[qid] = others
		if writer is not None:
			writer.writerow(aline)

                dict_germ_aln[qid] = my_alignment
		
		if my_alignment.sid not in dict_germ_count:
			dict_germ_count[my_alignment.sid] 	= 0
		dict_germ_count[my_alignment.sid] 		+= 1
		
		#total		+= 	1

	return dict_germ_aln, dict_j_seconds
 

def find_cdr3_borders(v_id,vgene,vlength,vstart,vend,jgene,jstart,j_start_on_read,jgaps,read_sequence):
	
	'''
	v_id = name of assigned V gene (eg IGHV1-2*02)
	vgene = germline sequence of assigned V gene
	vlength = length of QUERY sequence taken up by match (might be different from blast-reported length and/or vend-vstart+1 because of in-dels)
	vstart = position on germline V gene where match begins (hopefully = 1)
	vend = position on germline V gene where match ends
	jgene = germline sequence of assigned J gene
	jstart = position on germline J gene where match begins
	j_start_on_read = position on query (v-cut version, not full 454 read) where match with germline J begins
	jgaps = blast-reported number of gaps in J assignment
	read_sequence = V(D)J-trimmed sequence of the 454 read
	'''

	#print v_id,vgene,vlength,vstart,vend,jgene,jstart,j_start_on_read, v_align_len, read_sequence
	vMatches = []
	cys_pat = "TG[T|C|N]" #N is for a couple of shorter V's, like VH4-31
	if re.match("IGLV2-(11|23)", v_id):
		cys_pat = "TGCTGC" #special case
	if re.match("IGHV1-C",v_id):
		cys_pat = "TATGC"
	for m in re.finditer(cys_pat,vgene):
		vMatches.append(m)

	#last one **IN FRAME** is the cysteine we want! (matters for light chains)
	cdr3_start=-1
	vMatches.reverse()
	for cys in vMatches:
		if cys.start() % 3 == 0:
			cdr3_start = vlength - (vend - cys.start())# - (vstart - 1)
			break

	if cdr3_start > vlength:
		#uh oh, there is an indel in this read that caused blast to truncate the v gene alignment before reaching the cdr3
		#let's see if we can find the CxR/K motif to get the right place to start
		# (only works for heavy chain...)
		has_cxrk, c_start, c_end = has_pat(read_sequence[vlength:], pat_nuc_cxrk)
		if has_cxrk:
			cdr3_start = vlength + c_start
		else:
			#nothing more we can do right now
			pass

	jMotif = "TGGGG"
	if is_light:
		jMotif = "TT[C|T]GG"
	jMatch = re.search(jMotif,jgene)
	
	cdr3_end = vlength + j_start_on_read + (jMatch.start() - jstart) +3

	if jgaps > 0:
		#check for jMotif on read to correct for gaps
		for wgxg in re.finditer(jMotif, read_sequence):
                        #could reverse this as with cys above, but gap-cap takes care of other occurences
			if abs(wgxg.start() + 3 - cdr3_end) <= jgaps:
				#print "Found a TGGG motif at %d!" % wgxg.start()
				cdr3_end = wgxg.start() + 3
				break
	
	#print vlength, j_start_on_read,jMatch.start(), jstart
	#print cdr3_start, cdr3_end
	return cdr3_start, cdr3_end


def main():

	outfile = "%s/%s_VJtrim.fa" %(prj_tree.filtered, prj_name)
        handle = open(outfile, "w")

	v_only_out = "%s/%s_5trim.fa" %(prj_tree.filtered, prj_name)
        vhandle = open(v_only_out, "w")

	plus_bad_vj_out = "%s/%s_allJ.fa" %(prj_tree.filtered, prj_name)
        jhandle = open(plus_bad_vj_out, "w")

	cdr3_outfile = "%s/%s_cdr3.fa" %(prj_tree.filtered, prj_name)
        cdr3_handle = open(cdr3_outfile, "w")

	aa_outfile = "%s/%s_AAcdr3.fa" %(prj_tree.filtered, prj_name)
        aa_handle = open(aa_outfile, "w")

	plus_bad_outfile = "%s/%s_allCDR3.fa" %(prj_tree.filtered, prj_name)
        bad_handle = open(plus_bad_outfile, "w")

	print "curating 3' end..."
	#os.system("rm %s/%s_jgerm_tophit.txt" %(prj_tree.data, prj_name))

	total, found, noV, noJ, f_ind  = 0, 0, 0, 0, 1
	counts = {'good':0,'indel':0,'noCDR3':0,'stop':0,'out-of-frame':0}

	writer = csv.writer(open("%s/%s_jgerm_tophit.txt" %(prj_tree.data, prj_name), "w"), delimiter = sep)
	writer.writerow(PARSED_BLAST_HEADER)

	while os.path.isfile("%s/vcut_%06d.fasta" % (prj_tree.split, f_ind)):

		dict_germ_aln, dict_j_seconds =  get_top_hit("%s/vc%06d.txt" % (prj_tree.germ, f_ind), dict_j, writer)
		v_assignments, v_seconds = get_top_hit("%s/%s_%06d.txt" % (prj_tree.germ, prj_name, f_ind), dict_v)

		for entry in SeqIO.parse(open("%s/%s_%06d.fasta" % (prj_tree.split, prj_name, f_ind), "rU"), "fasta"):
			total += 1

			seq_id = str(int(entry.id))
			if not seq_id in v_assignments:
				noV+=1
			elif not seq_id in dict_germ_aln:
				noJ+=1
				myV = v_assignments[seq_id]
				if (myV.strand == '+'):
					entry.seq = entry.seq[ myV.qstart - 1 :  ]
				else:
					entry.seq = entry.seq[  : myV.qend ].reverse_complement()
				myVgenes = myV.sid
				if seq_id in v_seconds:
					myVgenes = myVgenes + "," + ",".join(v_seconds[seq_id])
				entry.description = "V_gene=%s status=noJ" % (myVgenes)
				vhandle.write(">%s %s\n" %(entry.id, entry.description))
				vhandle.write("%s\n" %entry.seq)
			else:
				found += 1

				myV = v_assignments[seq_id]
				myJ = dict_germ_aln[seq_id]
				status = "good"
				
				#if the 5' end is highly mutated, blast might leave it out. let's try to get it back
				if (myV.strand == "+" and myV.sstart>1 and myV.qstart>1):
					shift_back = min([myV.sstart, myV.qstart]) - 1
					#but start at an even codon
					shift_back -= (1 - ((myV.sstart - shift_back) % 3)) % 3
					addition = entry.seq[myV.qstart-shift_back : myV.qstart+2]
					if '*' not in addition.translate():
						myV.sstart    -= shift_back
						myV.qstart    -= shift_back
						myV.alignment += shift_back
				elif (myV.strand == "-" and myV.send>1 and myV.qend<len(entry.seq)):
					shift_back = min([myV.send-1, len(entry.seq)-myV.qend])
					#but start at an even codon
					shift_back -= (1 - ((myV.send - shift_back) % 3)) % 3
					addition = entry.seq[myV.qend-2 : myV.qend + shift_back].reverse_complement()
					if '*' not in addition.translate():
						myV.send      -= shift_back
						myV.qend      += shift_back
						myV.alignment += shift_back

				#check 3' end of J
				v_len   = myV.qend - (myV.qstart-1) #need to use qstart and qend instead of alignment to account for gaps
				#if (myJ.send < dict_j[myJ.sid].seq_len-1 and len(entry.seq) > v_len + myJ.qend):
				#	toAdd = min([dict_j[myJ.sid].seq_len-1-myJ.send,len(entry.seq)-v_len + myJ.qend])
				#	myJ.qend += toAdd
				if (myJ.send == dict_j[myJ.sid].seq_len):
					myJ.qend -= 1

				#get actual V(D)J sequence
				vdj_len = v_len + myJ.qend
				if (myV.strand == '+'):
					entry.seq = entry.seq[ myV.qstart - 1 : myV.qstart + vdj_len - 1 ]
				else:
					entry.seq = entry.seq[ myV.qend - vdj_len : myV.qend ].reverse_complement()

				#get CDR3 boundaries
				cdr3_start,cdr3_end = find_cdr3_borders(myV.sid,dict_v[myV.sid].seq, v_len, min(myV.sstart, myV.send), max(myV.sstart, myV.send), dict_j[myJ.sid].seq, myJ.sstart, myJ.qstart, myJ.gaps, entry.seq.tostring()) #min and max statments take care of switching possible minus strand hit
				cdr3_seq = entry.seq[ cdr3_start : cdr3_end ]
				#print cdr3_seq
				#break

				#check for continuous ORF
				v_frame = min([myV.sstart, myV.send]) % 3
				j_frame = 3 - ( (dict_j[myJ.sid].seq_len - myJ.sstart - 1) % 3 )
				#frame_shift = (myV.alignment + myJ.qstart - 1) % 3
				frame_shift = (v_len + myJ.qstart - 1) % 3
				if (v_frame + frame_shift) % 3 != j_frame % 3: #this can be wrong if there's a gap in the J alignment, but that doesn't seem to be frequent enough to be worth worrying about
					#frame +=1
					status = "out-of-frame"
					#continue

				#push the sequence into frame for translation, if need be
				five_prime_add = (v_frame-1) % 3
				entry.seq = 'N' * five_prime_add + entry.seq 

				#check for stop codons
				if '*' in entry.seq.translate():
					#stop+=1
					status="stop"
					#continue

				#even if V and J are in frame with no stop codons, the junction may still be out of frame due to gaps in the V alignment
				#theoretically, those gaps could let us do IMGT-style error correction, if we ask blast to print out the locations
				if len(cdr3_seq) % 3 != 0:
					#shift+=1
					status="indel"
					#continue

				#make sure cdr3 boundaries make sense
				if (cdr3_end<=cdr3_start or cdr3_end>vdj_len or cdr3_start<0):
					#print "No CDR3 for %s: start=%d, end=%d, read ends at %d"%(seq_id, cdr3_start, cdr3_end, vdj_len)
					#none+=1;
					status="noCDR3"
					#continue

				#add germline assignments to fasta description and write to disk
				myVgenes = myV.sid
				if seq_id in v_seconds:
					myVgenes = myVgenes + "," + ",".join(v_seconds[seq_id])
				myJgenes = myJ.sid
				if seq_id in dict_j_seconds:
					myJgenes = myJgenes + "," + ",".join(dict_j_seconds[seq_id])
				entry.description = "V_gene=%s J_gene=%s status=%s cdr3_nt_len=%d" % (myVgenes, myJgenes, status, len(cdr3_seq)-6)
				vhandle.write(">%s %s\n" %(entry.id, entry.description))
				vhandle.write("%s\n" %entry.seq)
				jhandle.write(">%s %s\n" %(entry.id, entry.description))
				jhandle.write("%s\n" %entry.seq)

				if status != "noCDR3":
					bad_handle.write(">%s %s\n" %(entry.id, entry.description))
					bad_handle.write("%s\n" % cdr3_seq)

				if status == "good":
			        	#may eventually want to add a length filter here, too
					#good += 1

					entry.description = "V_gene=%s J_gene=%s status=%s cdr3_len=%d" % (myVgenes, myJgenes, status, (len(cdr3_seq)/3)-2)
			        	#it's a high-quality sequence! Save to special file
					handle.write(">%s %s\n" %(entry.id, entry.description))
					handle.write("%s\n" %entry.seq)
					cdr3_handle.write(">%s %s\n" %(entry.id, entry.description))
					cdr3_handle.write("%s\n" % cdr3_seq)
					aa_handle.write(">%s %s\n" %(entry.id, entry.description))
					aa_handle.write("%s\n" % cdr3_seq.translate())

				counts[status] += 1

		print "%d done, found %d; %d good..." %(total, found, counts['good'])
		f_ind += 1

	#print "%d done, found %d; %d good." %(total, found, counts['good'])
	print "\tNo V assignment: %d\n\tNo J assignment: %d\n\tNo CDR3 found: %d\n\tCDR3 Indels: %d\n\tStop codons (anywhere): %d\n\tOut-of-frame junction (despite a good CDR3 and no stop codons): %d\n"%(noV,noJ,counts['noCDR3'],counts['indel'],counts['stop'],counts['out-of-frame'])
	handle.close()

	#now get original raw reads and write an overall stat file
	raw=0
        fastas = glob.glob("%s/*.fa"  %prj_tree.original) + glob.glob("%s/*.fasta" %prj_tree.original) + glob.glob("%s/*.fna" %prj_tree.original)
        for fasta_file in fastas:
		for entry in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
			raw += 1


	write_dict2file(dict_germ_count, total, ["subject", "count", "percent"])
	stats_file = "%s/%s_read_stats.tab" %(prj_tree.data,prj_name)
	stats_writer = csv.writer(open(stats_file, "w"), delimiter = sep)
	stats_writer.writerow(["Raw","Correct_Length","V_assigned","J_assigned","CDR3_assigned","ORF"])
	stats_writer.writerow([raw,total,total-noV,total-noV-noJ,total-noV-noJ-counts['noCDR3'],counts['good']])
	stats_writer.writerow([])
	stats_writer.writerow(["CDR3 Indels:",counts['indel']])
	stats_writer.writerow(["Stop codons (anywhere):",counts['stop']])
	stats_writer.writerow(["Out-of-frame junctions (despite a good CDR3 and no stop codons):",counts['out-of-frame']])



if __name__ == '__main__':
	
	# get parameters from input
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, l="is_light")
	is_light = getParas(dict_args, "is_light")	
	
	prj_tree 		= ProjectFolders(os.getcwd())
	prj_name 		= fullpath2last_folder(prj_tree.home)

	v_db,j_db      =  dict_germ_db[is_light], dict_jgerm_db[is_light]
	dict_v         =  load_fastas(v_db)
	dict_j         =  load_fastas(j_db)
	dict_strand_count 		= {"+" : 0.0, "-" : 0.0}
	dict_germ_count			= dict()
	


	main()

