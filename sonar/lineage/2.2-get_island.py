#!/usr/bin/env python

"""
2.2-get_island.py

This script selects sequences based on their positions in an identity-divergence
      plot and copies them to a unique file.

Usage: 2.2-get_island.py -n native [-imin min_identity -imax max_indentity
                                    -dmin min_divergence -dmax max_divergence
				    -seq seqs.fa -div id-div.tab]

    Invoke with -h or --help to print this documentation.

    native	Name of the known mAb used as the identity referent.
    imin	Minimum percent sequence identity to the identity referent.
                   Default = 85.
    imax	Maximum percent sequence identity to the identity referent.
                   Included for completeness. Default = 100.
    dmin	Minimum percent germline divergence. Default = 0.
    dmax	Maximum percent germline divergence. Default = 40.
    seq		Custom input file with sequences to search. By default, uses
                   pipeline output with "goodVJ_unique" nucleotide sequences.
    div		If a custom sequence file is used, please specify the location
                   of the corresponding output from 2.1-calculate_id-div.pl.

Created by Chaim A Schramm, 2012-10-04.
Edited and commented for publication by Chaim A Schramm on 2015-04-20.

Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys

try:
        from sonar.lineage import *
except ImportError:
        find_SONAR = sys.argv[0].split("sonar/annotate")
        sys.path.append(find_SONAR[0])
        from sonar.lineage import *

global inFile, native


def main():
	
	if not os.path.isfile(divFile):
		sys.exit("Please run 2.1-calulate_id-div.pl before running this script!")

	reader 	= csv.reader(open(divFile, "rU"), delimiter = sep)
	natives = reader.next()[2:]

	try:
		pos = natives.index(native)
	except:
		sys.exit("Can't find desired mAb in %s (options are %s)." % (divFile, ", ".join(natives)))

	# a dictionary for keeping track of those that survive the filter
	island = []
		
	bad = 0
	for row in reader:
		if len(row) - 2 < len(natives):
			pass;
		elif row[1] == "NA" or row[2+pos] == "NA":
			bad+=1
		else:
			read_id, divergence, identity = row[0], float(row[1]), float(row[2+pos])
			if (divergence >= min_div and divergence <= max_div and identity >= min_iden and identity <= max_iden):
				island.append(read_id)


	#error checking
	if len(island) == 0:
		sys.exit( "No reads were in found in the specificed island (referent = %s; id: %d-%d; div: %d-%d).\n Please check boundaries and try again." % 
			  ( native, min_iden, max_iden, min_div, max_div ) )

	#now get the sequences
	print "Found %d sequences in the island" % len(island)
	island_reads = load_fastas_in_list(inFile, island)

	#and write
	outfile = "%s/%s_%s_id%d-%d_div%d-%d.fasta" % (prj_tree.nt, prj_name, native, min_iden, max_iden, min_div, max_div)
	handle = open(outfile, "w")
	SeqIO.write(island_reads.values(), handle, "fasta")
	handle.close()



if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	if q("-seq") and not q("-div"):
		print "Please specify location of output from 2.1-calculate_id-div.pl for custom sequence file!"
		sys.exit(1)

	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	= fullpath2last_folder(prj_tree.home)
	
	dict_args = processParas(sys.argv, n="native",imin="min_iden",imax="max_iden",dmin="min_div",dmax="max_div", seq="inFile", div="divFile")
	defaults = dict( native="",min_iden=85,max_iden=100,min_div=0,max_div=40,inFile="%s/%s_goodVJ.fa"%(prj_tree.nt, prj_name),divFile="%s/%s_goodVJ_unique_id-div.tab" % (prj_tree.tables, prj_name) )
	native,min_iden,max_iden,min_div,max_div,inFile,divFile = getParasWithDefaults(dict_args, defaults, "native","min_iden","max_iden","min_div","max_div","inFile", "divFile")
	
	if native == "":
		print "Please specify native mAb to use as identity referent.\n\n"
		print __doc__
		sys.exit(1)

	main()

