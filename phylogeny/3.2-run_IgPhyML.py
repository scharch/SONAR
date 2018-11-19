#!/usr/bin/env python3

"""
3.2-runIgPhyML.py

This script calls IgPhyML to generate a maximum likelihood tree for a set
      of sequences, rooted on the germline V gene sequence. It then also
      calls the IgPhyl Utilities script for ancestor recontruction. For
      optimal results, sequences should be aligned manually and specified
      with the -i parameter. However, the script can also use MUSCLE to
      create an automated alignment. For automatic alignments, please use 
      the -v parameter to specify the assigned V gene, which will be added
      to the alignment and used for rooting the tree. 
Tree file goes in output/<prj_name>_igphyml.tree; inferred sequences go in
      output/sequences/(nucleotide|amino_acid)/<prj_name>_inferredAncestors.fa;
      other IgPhyML output goes in output/logs/<prj_name>_igphyml_stats.txt.

For more information about IgPhyML, please see
      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5419485/
and
      https://github.com/kbhoehn/IgPhyML


Usage: 3.2-runIgPhyML.py -i input [ --format fasta --root IGHV1-2*02 --quick -f ]
       3.2-runIgPhyML.py -v IGHV3-30*18 [ --locus <H|K|L> | --lib path/to/library.fa ] [ --seqs input.fa --natives natives.fa --quick -f ]

Options:
    -i input                   Manual alignment (in PHYLIP format) of the sequences to be 
                                  analayzed, with known antibody seqeunces and germline 
                                  (or other outgroup) sequence included. This is the 
                                  preferred option for running this program and should be
                                  used especially for inferring ancestral sequences.
                                  If --natives and -v (see below) are supplied instead, 
                                  sequences will be aligned automagically with MUSCLE.
    --format fasta             Format of the alignment, eg fasta or phylip. [default: fasta]
    --root IGHV1-2*02          Name of the sequence to use as outgroup/root of the tree.
                                  Typically detected automatically if it's a germline V
                                  gene, but may need to be specified if eg, your using an
                                  estimated naive sequence as the root.
    --quick                    Flag to speed up calculation by turning off SPR moves when 
                                  finding initial tree and skipping topology refinement when 
                                  recalulating with the HLP17 model. [default: False]
    -v IGHV3-30*18             Assigned germline V gene of known antibodes, for use in 
                                  rooting the trees. Include allele designation.
    --locus <H|K|L>            Specify use of V heavy/kappa/lambda germlines libraries,
                                  respectively. Mutually exclusive with --lib. [default: H]
    --lib path/to/library.fa   Optional custom germline library (eg for Rhesus or Mouse).
    --seqs input.fa            A fasta file containing the sequences from which the tree is
                                  to be built. [default: output/sequences/nucleotide/<project>-collected.fa]
    --natives natives.fa       A fasta file containing known sequences to be included in
                                  the tree.
    -f		               Force a restart of the analysis, even if there are files from
                                  a previous run in the working directory.

Created by Chaim A Schramm 2015-07-09 as 3.2-run_DNAML.py.
Added options for flexibility and more informative error messages by CAS 2016-08-19.
Edited to use Py3 and DocOpt by CAS 2018-08-29.
Changed tree-building engine to IgPhyML and renamed by CAS 2018-10-22.
Added gap handling (ported from 3.3) by CAS 2018-10-23.

Copyright (c) 2011-2018 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from collections import defaultdict, OrderedDict
from docopt import docopt
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

try:
	from sonar.phylogeny import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/phylogeny")
	sys.path.append(find_SONAR[0])
	from sonar.phylogeny import *



def deleteGap( start, end, name, tree, gaps ):

	for i,g in gaps[name]:
		if g['start']<= start and g['end']>=end:
			if g['value'] > .5:
				return #wasn't provisional this far down

			#check boundaries to see if we need to split
			if g['end']>end:
				gaps[name].insert(i,{'start':end, 'end':g['end'], 'value':g['value']})
			if g['start']<start:
				gaps[name].insert(i,{'start':g['start'], 'end':start, 'value':g['value']})

			#delete
			gaps[name].remove(g)
			break

	#iterate down the tree
	for child in tree[name]['children']:
		deleteGap(start,stop,child,tree,gaps)



def assignGaps( name, tree, gaps ):

	provisional = []
	for child in tree[name]['children']:
		assignGaps(child, tree, gaps)
		for g in gaps[child]:
			provisional.append(g.copy())

	sortProv = sorted( provisional, key=lambda k: (k['start'], k['end']) )
	ind = 0
	while ind < len(sortProv):
		if ind+1 < len(sortProv):
			if sortProv[ind]['end']>sortProv[ind+1]['start']:
				#these gaps overlap
				s1, e1, v1 = sortProv[ind]['start'], sortProv[ind]['end'], sortProv[ind]['value']
				s2, e2, v2 = sortProv[ind+1]['start'], sortProv[ind+1]['end'], sortProv[ind]['value']

				if s2>s1:
					gaps[name].append( { 'start':s1, 'end':s2, 'value':(v1+v2)/2 } )
				gaps[name].append( { 'start':s2, 'end':e1, 'value':(v1+v2)/2 } )
				if e2>e1:
					gaps[name].append( { 'start':e1, 'end':e2, 'value':(v1+v2)/2 } )
				
				ind += 2 #skip the overlapping one we just processed

			else:
				#this one stands by itself
				gaps[name].append( sortProv[ind].copy() )
				ind += 1

		else:
			#last gap, stands by itself
			gaps[name].append( sortProv[ind].copy() )
			ind += 1

	#now that we've gone through all of the descendant gaps, remove ones that aren't sufficiently supported
	for gg in gaps[name]:
		if gg['value']<0.5:
			deleteGap(gg['start'], gg['end'], name, tree, gaps)



def getFinalSeqs( seqDict, gapList, trans=False ):
	for seqID in seqDict:
		myseq = seqDict[seqID]
		for g in gapList[seqID]:
			myseq.seq= myseq.seq[0:g['start']] + "-" * (g['end']-g['start']) + myseq.seq[g['end']:]
		if trans:
			myseq.seq = myseq.seq.translate(table=GAPPED_CODON_TABLE)
		yield myseq



def main():

	oldFiles = glob.glob("%s/infile"%prj_tree.phylo) + glob.glob("%s/%s_igphyml.tree"%(prj_tree.out,prj_name)) + glob.glob("%s/%s_igphyml_stats.txt"%(prj_tree.logs,prj_name))
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
		run_muscle	      = MuscleCommandline( cmd=muscle, input="%s/%s_to_align.fa" % (prj_tree.phylo, prj_name), out="%s/%s_aligned.afa" % (prj_tree.phylo, prj_name) )
		run_muscle.maxiters   = 2
		run_muscle.diags      = True
		run_muscle.gapopen    = -5000.0 #code requires a float
		print( run_muscle )
		run_muscle()

		#this is probably bad form
		arguments['-i'] = "%s/%s_aligned.afa" % (prj_tree.phylo, prj_name)


	#open the alignment to rename everything and find germline sequence
	with open(arguments['-i'], "rU") as handle:
		try:
			aln = AlignIO.read(handle, arguments['--format'])
		except:
			sys.exit("Couldn't read alignment: is %s the correct format?"%arguments['--format'])
				
	align_len = aln.get_alignment_length()
	
	#kill the fasta def line and any usearch/vsearch annotations to avoid formatting foul-ups
	germ_id	  = ""
	foundRoot = False
	gaps	  = defaultdict( list )
	for seq in aln:
		seq.id = re.sub("[;:].*", "", seq.id)
		seq.description = ""
		if re.search("(IG|VH|VK|VL|HV|KV|LV)", seq.id, re.I) is not None:
			germ_id = seq.id
			
		if arguments['--root'] is not None and seq.id == arguments['--root']:
			foundRoot = True
			
		for g in re.finditer("-+", str(seq.seq)):
			#save gap. value is a field to help me determine what's real in assignGaps
			gaps[ seq.id ].append( {'start':g.start(), 'end':g.end(), 'value':1} )
			  
	if arguments['--root'] is not None:
		germ_id = arguments['--root']
		if not foundRoot:
			sys.exit("Couldn't find specified root sequence %s in input file"%arguments['--root'])
	elif germ_id=="":
		sys.exit("Couldn't find a germline gene in the alignment, please use the --root option and try again.")
		
	with open("%s/infile" % prj_tree.phylo, "w") as output:
		AlignIO.write(aln, output, "fasta")


	#now call IgPhyML
	#fast initial tree
	opts = ["-s", "SPR"]
	if arguments['--quick']:
		opts=[]
	#set an environmental variable so that IgPhyML can find its libraries
	os.environ.update( { 'IGPHYML_PATH' : '%s/third-party/src/motifs'%SCRIPT_FOLDER } )
	s = subprocess.Popen([igphyml, "-i", "%s/infile" % prj_tree.phylo,
			      "-m", "GY", "-w", "MO", "-t", "e", "--run_id", "gy94"] + opts,
			     universal_newlines=True, stderr=subprocess.PIPE)
	o,e = s.communicate()
	
	if re.search("error while loading shared libraries", str(e)):
		#Some libraries needed for optimized execution are missing
		#  Try again with a version compiled without optimizations
		s = subprocess.Popen([igphyml_slow, "-i", "%s/infile" % prj_tree.phylo,
				      "-m", "GY", "-w", "MO", "-t", "e", "--run_id", "gy94"] + opts,
				     universal_newlines=True, stderr=subprocess.PIPE)
		o,e = s.communicate()

	if e != "" or s.returncode != 0:
		sys.exit( "Error running '%s':\n%sExit code %d" % (" ".join(s.args),e,s.returncode) )
		

	#Refine tree with AID-specific hotpsot motifs
	opts = ['-o', 'tlr']
	if arguments['--quick']:
		opts = ['-o', 'lr']
	s = subprocess.Popen([igphyml, "-i", "%s/infile" % prj_tree.phylo,
			      "-m", "HLP17", "--root", germ_id,
			      "-u", "%s/infile_igphyml_tree.txt_gy94" % prj_tree.phylo,
			      "--motifs", "FCH", "--run_id", "hlp17",
			      "--ambigfile", "%s/ambigfile.txt" % prj_tree.phylo] + opts,
			     universal_newlines=True, stderr=subprocess.PIPE)
	o,e = s.communicate()

	if re.search("error while loading shared libraries", str(e)):
		s = subprocess.Popen([igphyml_slow, "-i", "%s/infile" % prj_tree.phylo,
				  "-m", "HLP17", "--root", germ_id,
				  "-u", "%s/infile_igphyml_tree.txt_gy94" % prj_tree.phylo,
				  "--motifs", "FCH", "--run_id", "hlp17",
				  "--ambigfile", "%s/ambigfile.txt" % prj_tree.phylo] + opts,
				 universal_newlines=True, stderr=subprocess.PIPE)
		o,e = s.communicate()

	if e != "" or s.returncode != 0:
		sys.exit( "Error running '%s':\n%sExit code %d" % (" ".join(s.args),e,s.returncode) )
		

	#now need to set up a config file for ancestor reconstruction
	with open("%s/ar.config"%prj_tree.phylo, "w") as handle:
		handle.write( "length\t%d\n" % (align_len/3) )
		handle.write( "rooted\t1\noutdir\t%s\n" % prj_tree.phylo )
		handle.write( "seqfile\t%s/infile\n" % prj_tree.phylo )
		handle.write( "rootid\t%s\n" % germ_id )
		handle.write( "igphyml\t%s/%s\n" % (SCRIPT_FOLDER, "third-party") )
		handle.write( "stats\t%s/infile_igphyml_stats.txt_hlp17\n" % prj_tree.phylo )
		handle.write( "tree\t%s/infile_igphyml_tree.txt_hlp17\n" % prj_tree.phylo )
		handle.write( "ambigfile\t%s/ambigfile.txt\n" % prj_tree.phylo )
		handle.write( "stem\t%s\n" % prj_name )

	s = subprocess.Popen( ["perl", "-I", "%s/third-party"%SCRIPT_FOLDER, reconstruct, "%s/ar.config"%prj_tree.phylo], universal_newlines=True, stderr=subprocess.PIPE )
	o,e = s.communicate()
	if e != "" or s.returncode != 0:
		sys.exit( "Error running '%s':\n%sExit code %d" % (" ".join(s.args),e,s.returncode) )


	if len(gaps)>0:
		#fix ancestor inference by putting gaps back in
		#start by reading in inferred sequences and reconstructing the tree
		tree	= dict()
		stack	= list()
		seqDict = OrderedDict()
		with open("%s/%s.MLcodons.fa"%(prj_tree.phylo,prj_name), "rU") as infer:
			for seq in SeqIO.parse(infer, "fasta"):
				name = seq.id.split(";")[1]
				seqDict[name] = seq
				if "," in name:
					tree[name] = { 'id':name, 'children':stack[-2:] }
					tree[stack.pop()]['parent'] = name
					tree[stack.pop()]['parent'] = name
					stack.append( name )
				else:
					tree[name] = { 'id':name, 'children':[] }
					stack.append( name )

		#now iterate down tree to propogate gaps
		assignGaps(stack[0], tree, gaps)

		#do output
		with open( "%s/%s_inferredAncestors.fa"%(prj_tree.nt,  prj_name), "w" ) as handle:
			SeqIO.write( getFinalSeqs(seqDict,gaps), handle, "fasta" )
		with open( "%s/%s_inferredAncestors.fa"%(prj_tree.aa,  prj_name), "w" ) as handle:
			SeqIO.write( getFinalSeqs(seqDict,gaps, trans=True), handle, "fasta" )
		
	else:
		os.rename( "%s/%s.MLcodons.fa"%(prj_tree.phylo,prj_name), "%s/%s_inferredAncestors.fa"%(prj_tree.nt,prj_name) )
		os.rename( "%s/%s.MLaas.fa"   %(prj_tree.phylo,prj_name), "%s/%s_inferredAncestors.fa"%(prj_tree.aa,prj_name) )

		
	#move non-seqeunce outputs to logical places
	os.rename( "%s/infile_igphyml_stats.txt_hlp17"%prj_tree.phylo,		  "%s/%s_igphyml_stats.txt"   %(prj_tree.logs,prj_name) )
	os.rename( "%s/infile_igphyml_tree.txt_hlp17" %prj_tree.phylo,		  "%s/%s_igphyml.tree"	      %(prj_tree.out, prj_name) )


if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)
	natives	 = {} #avoid errors
	
	
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


	if arguments['--root'] is not None:
		if re.search("[;:]",arguments['--root']):
			arguments['--root'] = re.sub("[;:].*","",arguments['--root'])
			sys.stderr.write("Found forbidden character (';' or ':') in root id, will use %s instead"%arguments['--root'])
			
	#check directories to avoid errors
	os.makedirs(prj_tree.phylo, exist_ok=True)
	os.makedirs(prj_tree.aa, exist_ok=True)
	os.makedirs(prj_tree.nt, exist_ok=True)
	os.makedirs(prj_tree.logs, exist_ok=True)
		
	#log command line
	logCmdLine(sys.argv)
	
	
	main()

