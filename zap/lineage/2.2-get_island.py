#!/usr/bin/python

"""
2.2-get_island.py

This script selects sequences based on their positions in an identity-divergence
      plot and copies them to a unique file.

Usage: 2.2-get_island.py -n native [-imin min_identity -imax max_indentity
                                    -dmin min_divergence -dmax max_divergence
				    -i seqs.fa]

    Invoke with -h or --help to print this documentation.

    native	Name of the known mAb used as the identity referent.
    imin	Minimum percent sequence identity to the identity referent.
                   Default = 85.
    imax	Maximum percent sequence identity to the identity referent.
                   Included for completeness. Default = 100.
    dmin	Minimum percent germline divergence. Default = 0.
    dmax	Maximum percent germline divergence. Default = 40.
    i		Custom input file with sequences to search. By default, uses
                   pipeline output with "goodVJ" nucleotide sequences.

Created by Chaim A Schramm, 2012-10-04.
Edited and commented for publication by Chaim A Schramm on 2015-04-20.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

from zap import *

global inFile, native


def main():
	
	if not os.path.isfile("%s/%s_id-div.tab" % (prj_tree.tables, prj_name)):
		sys.exit("Please run 2.1-calulate_id-div.pl before running this script!")

	reader 	= csv.reader(open("%s/%s_id-div.tab" % (prj_tree.tables, prj_name), "rU"), delimiter = sep)
	natives = reader.next()[2:]

	try:
		pos = natives.index(native)
	except:
		sys.exit("Can't find desired mAb in %s/%s_id-div.tab (options are %s)." % (prj_tree.tables, prj_name, ", ".join(natives)))

	# a dictionary for keeping track of those that survive the filter
	island = []
		
	for row in reader:
		if len(row) - 2 < len(natives):
			pass;
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

	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	= fullpath2last_folder(prj_tree.home)
	
	dict_args = processParas(sys.argv, n="native",imin="min_iden",imax="max_iden",dmin="min_div",dmax="max_div", i="inFile")
	defaults = dict( native="",min_iden=85,max_iden=100,min_div=0,max_div=40,inFile="%s/%s_goodVJ.fa"%(prj_tree.nt, prj_name) )
	native,min_iden,max_iden,min_div,max_div,inFile = getParasWithDefaults(dict_args, defaults, "native","min_iden","max_iden","min_div","max_div","inFile")
	
	if native == "":
		print "Please specify native mAb to use as identity referent.\n\n"
		print __doc__
		sys.exit(1)

	main()

