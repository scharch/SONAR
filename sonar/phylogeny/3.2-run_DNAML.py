#!/usr/bin/env python3

"""
3.2-runDNAML.py

This script calls DNAML to generate a maximum likelihood tree for a set
      of sequences, outgroup-rooted on the germline V gene sequence. For
      optimal results, sequences should be aligned manually (use the 
      DNAML/PHYLIP format) and specified with the -i parameter. However,
      the script can also use MUSCLE to create an automated alignment to
      be passed directly to DNAML. For automatic alignments, please use 
      the -v parameter to specify the assigned V gene, which will be added
      to the alignment and used for rooting the tree.


Usage: 3.2-runDNAML.py -i input.phy [ --outtree out.tree --outfile dnaml.out -j 5 -g -f ]
       3.2-runDNAML.py -v IGHV3-30*18 [ --locus <H|K|L> | --lib path/to/library.fa ] [ --seqs input.fa --natives natives.fa ] [ --outtree out.tree --outfile dnaml.out -j 5 -g -f ]

Options:
    -i input.phy               Manual alignment (in PHYLIP format) of the sequences to be 
                                  analayzed, with known antibody seqeunces and germline 
                                  (or other outgroup) sequence included. This is the 
                                  preferred option for running this program and should be
                                  used especially for inferring ancestral sequences.
                                  If --natives and -v (see below) are supplied instead, 
                                  sequences will be aligned automagically with MUSCLE.
    -v IGHV3-30*18             Assigned germline V gene of known antibodes, for use in 
                                  rooting the trees. Include allele designation.
    --locus <H|K|L>            Specify use of V heavy/kappa/lambda germlines libraries,
                                  respectively. Mutually exclusive with --lib. [default: H]
    --lib path/to/library.fa   Optional custom germline library (eg for Rhesus or Mouse).
    --seqs input.fa            A fasta file containing the sequences from which the tree is
                                  to be built. [default: output/sequences/nucleotide/<project>-collected.fa]
    --natives natives.fa       A fasta file containing known sequences to be included in
                                  the tree.
    --outtree out.tree         Where to save the output tree. [default: output/<project>.tree]
    --outfile dnaml.out        Where to save DNAML output (text tree and inferred ancestors)
                                  [default: output/logs/<project>.dnaml.out]
    -j 5                       Number of times to "jumble" the order of the input sequences
                                  and rebuild the tree [default: 5]
    -g                         Flag to have DNAML do global rearrangments. Can make things
                                  v e r y  s l o w for large alignments. [default: False]
    -f		               Force a restart of the analysis, even if there are files from
                                  a previous run in the working directory.

Created by Chaim A Schramm 2015-07-09.
Added options for flexibility and more informative error messages by CAS 2016-08-19.
Edited to use Py3 and DocOpt by CAS 2018-08-29.

Copyright (c) 2011-2018 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

try:
	from sonar.phylogeny import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/phylogeny")
	sys.path.append(find_SONAR[0])
	from sonar.phylogeny import *



def revertName(match):
	return lookup[ int(match.group(0)) - 1 ] #ids are 1-indexed, list is 0-indexed



def main():

	global lookup

	oldFiles = glob.glob("%s/infile"%prj_tree.phylo) + glob.glob("%s/outtree"%prj_tree.phylo) + glob.glob("%s/outfile"%prj_tree.phylo)
	if len(oldFiles) > 0:
		if arguments['-f']:
			for f in oldFiles:
				os.remove(f)
		else:
			sys.exit("Old files exist! Please use the -f flag to force overwrite.")
	

	if arguments['-v'] is not None:

		#do alignment

		#first create a working file to align and add the germline and natives
		shutil.copyfile(arguments['--seqs'], "%s/%s_to_align.fa"%(prj_tree.phylo, prj_name))
		handle = open( "%s/%s_to_align.fa"%(prj_tree.phylo, prj_name), "a" )
		handle.write( ">%s\n%s\n" % (germ_seq.id, germ_seq.seq) )
		for n in natives.values():
			handle.write( ">%s\n%s\n" % (n.id, n.seq) )
		handle.close()

		#now run muscle
		run_muscle	      = MuscleCommandline( input="%s/%s_to_align.fa" % (prj_tree.phylo, prj_name), out="%s/%s_aligned.afa" % (prj_tree.phylo, prj_name) )
		run_muscle.maxiters   = 2
		run_muscle.diags      = True
		run_muscle.gapopen    = -5000.0 #code requires a float
		print( run_muscle )
		run_muscle()

		#this is probably bad form
		arguments['-i'] = "%s/%s_aligned.afa" % (prj_tree.phylo, prj_name)


	#open the alignment to rename everything and find germline sequence
	#rename is to avoid possible errors with DNAML from sequence ids that are too long
	germ_pos = 1
	with open(arguments['-i'], "rU") as handle:
		if arguments['-v'] is not None:
			aln = AlignIO.read(handle, "fasta")
		else: 
			try:
				aln = AlignIO.read(handle, "phylip-relaxed")
			except:
				sys.exit("Please make sure custom input is aligned and in PHYLIP format...")

				
	#change all ids to 10 digit numbers to avoid formatting foul-ups
	lookup = []
	for seq in aln:
		lookup.append( seq.id )
		if re.search("(IG|VH|VK|VL|HV|KV|LV)", seq.id) is not None:
			germ_pos = len( lookup )
		seq.id = "%010d" % len( lookup )


	with open("%s/infile" % prj_tree.phylo, "w") as output:
		AlignIO.write(aln, output, "phylip")


	#now generate script for DNAML
	# J is "jumble" followed by random seed and number of times to repeat
	# G is to do global rearrangements
	# O is outgroup root, followed by position of the germline in the alignment
	# 5 tells DNAML to do the ancestor inference
	# Y starts the run
	with open("%s/dnaml.in"%prj_tree.phylo, "w") as handle:
		seed = random.randint(0,1e10) * 2 + 1 #seed must be odd
		handle.write( "J\n%d\n%d\n" % (seed, arguments['-j']) )
		if (arguments['-g']):
			handle.write("G\n")
		handle.write( "O\n%d\n5\nY\n" % germ_pos )

		
	# change to work directory so DNAML finds "infile" and puts the output where we expect
	origWD = os.getcwd()
	os.chdir(prj_tree.phylo)
	with open("dnaml.in", "rU") as pipe:
		subprocess.call([dnaml], stdin=pipe)
	os.chdir(origWD)

	#revert names in tree
	with open("%s/outtree"%prj_tree.phylo, "rU") as intree:
		mytree = intree.read()
	fixedtree = re.sub("\d{10}", revertName, mytree)
	with open(arguments['--outtree'], "w") as outtree:
		outtree.write(fixedtree)

	#revert names in out file
	with open("%s/outfile"%prj_tree.phylo, "rU") as instuff:
		mystuff = instuff.read()
	fixedstuff = re.sub("\d{10}", revertName, mystuff)
	with open(arguments['--outfile'], "w") as outstuff:
		outstuff.write(fixedstuff)
	
	
	print( "\nOutput in %s and %s\n" % (arguments['--outtree'], arguments['--outfile']) )


	

if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)
	natives	 = {} #avoid errors
	arguments['-j'] = int(arguments['-j'])
	
	
	#first decide if program was called with -i or -v
	if arguments['-i'] is not None:
		if not os.path.isfile(arguments['-i']):
			sys.exit("Can't find input file %s" % argument['-i'])
	else:
		#find the right germline library
		if arguments['--lib'] is not None:
			if not os.path.isfile(arguments['-lib']):
				sys.exit("Can't find germline file %s" % argument['-lib'])
		else:
			if arguments['--locus'] in dict_vgerm_db.keys():
				arguments['--lib'] = dict_vgerm_db[ arguments['--locus'] ]
			else:
				sys.exit("Error: --locus must be one of H, K, or L")

		germ_dict = load_fastas(arguments['--lib'])

		#now get the V gene
		if not arguments['-v'] in germ_dict:
			sys.exit( "Specified germline gene (%s) is not present in the %s library!\n" % (arguments['-v'], arguments['--lib']) )
		germ_seq = germ_dict[ arguments['-v'] ]

		#get sequences to align
		arguments['--seqs'] = re.sub( "<project>", prj_name, arguments['--seqs'] )
		if not os.path.isfile(arguments['--seqs']):
			sys.exit( "Can't find sequence file %s to build tree from." % arguments['--seqs'] )
		
		#now load mAb sequences, if provided	    
		if arguments['--natives'] is not None:
			natives = load_fastas( arguments['--natives'] )
		else:
			print( "No native sequences specified; tree will only include NGS sequences." )

	#check output paths
	arguments['--outtree'] = re.sub( "<project>", prj_name, arguments['--outtree'] )
	arguments['--outfile'] = re.sub( "<project>", prj_name, arguments['--outfile'] )

	#check directories to avoid errors
	os.makedirs(prj_tree.phylo, exist_ok=True)
	os.makedirs(prj_tree.nt, exist_ok=True)
	os.makedirs(prj_tree.logs, exist_ok=True)
		
	#log command line
	logCmdLine(sys.argv)
	
	
	main()

