
from .. import *


def parse_hits(hits, handle=None, maxQEnd=99999, minQStart=-1, strand=None):

	best_alignment = MyAlignment( hits[0].split("\t") )
	best_hit = hits[0]
	others = []
	second_match = ""

	for h in hits[1:]:

		if h=="":
			continue #prevent errors from empty lines at the end

		my_alignment = MyAlignment( h.split("\t") )

		#CAS 2020-01-02
		if strand is not None and not my_alignment.strand == strand:
			continue

		#skips D genes that matched 5' J
		if my_alignment.qend > maxQEnd:
			old_id=""
			continue

		#skips C genes that matched 3' J (not sure if necessary)
		if my_alignment.qstart < minQStart:
			old_id=""
			continue

		#added 20150107 by CAS
		'''
		need three conditions:
		1. hit is on same gene
		2. hit is on same strand
		3. hits are non-overlapping
		'''
		if my_alignment.sid == best_alignment.sid and my_alignment.strand == best_alignment.strand and (max(my_alignment.sstart, my_alignment.send)<min(best_alignment.sstart,best_alignment.send) or min(my_alignment.sstart, my_alignment.send)>max(best_alignment.sstart,best_alignment.send)):
			second_match = h
			#change boundaries of alignment on both query and hit (to get J properly)
			if my_alignment.send < best_alignment.sstart:
				best_alignment.sstart = my_alignment.sstart
				if best_alignment.strand == "plus":
					best_alignment.qstart = my_alignment.qstart
				else:
					best_alignment.qend = my_alignment.qend
			else:
				best_alignment.send = my_alignment.send
				if best_alignment.strand == "plus":
					best_alignment.qend = my_alignment.qend
				else:
					best_alignment.qstart = my_alignment.qstart

		elif re.match("IG[HKL]J", my_alignment.sid) and my_alignment.score>=40 and my_alignment.qstart<best_alignment.qstart:
			#a bit of kludge for double J matches. 
			#Usually these are bad amplicons (or bad assemblies from single cell data)
			#My assumption is that the one closer to the V is the more reliable one
			best_alignment = my_alignment
			best_hit = h
				
		elif my_alignment.score >= best_alignment.score - 3 and my_alignment.sid.split("*")[0] != best_alignment.sid.split("*")[0] and not any( my_alignment.sid.split("*")[0] == x.split("*")[0] for x in others ):
			others.append(my_alignment.sid)

	#parsed all hits, do output and return
	if handle is not None:
		handle.write(best_hit+"\n")
		if second_match != "":
			handle.write(second_match+"\n")

	best_alignment.others = ",".join(others)
	return best_alignment
