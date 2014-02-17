#!/usr/bin/env python
# encoding: utf-8
"""
clustal-consensus.py

Created by Zhenhai Zhang on 2013-02-28.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from mytools import *
import dendropy
nucleotides = ["A", "C", "G", "T", "N", "-"]
def get_score(seq, d):
	score = 0.0
	for nuc in nucleotides:
		is_nuc = [x == nuc for x in seq]
		wm = list(d[nuc])
		zs = zip(is_nuc, wm)
		nscores = [x * y for x, y in zs]
		score += sum(nscores)
	return score
		


def main():
	infiles = sorted(glob.glob("*.aln"))
	handle = open("G01-99-97-rep-consensus.fasta", "w")
	for infile in infiles:
		align = AlignIO.read(infile, "clustal")
		total, dict_ab, aln_len = len(align), dict(), align.get_alignment_length()
		for i in range(total):
			dict_ab[align[i].id] = align[i].seq.tostring()
			

		dict_nuc = dict()
		for nuc in nucleotides:
			dict_nuc[nuc] = zeros(aln_len)
			for k, v in dict_ab.items():
				is_nuc = array([x == nuc for x in v])
				dict_nuc[nuc] += is_nuc / float(total)
		max_score, max_id = 0.0, ""
		for k, v in dict_ab.items():
			score = get_score(v, dict_nuc)
			if score > max_score:
				max_score = score
				max_id = k
		SeqIO.write(SeqRecord(Seq.Seq("".join(dict_ab[max_id].split("-"))), id=max_id, description=max_id), handle, "fasta")
		
		
		


if __name__ == '__main__':
	

	main()

