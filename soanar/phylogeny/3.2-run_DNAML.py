#!/usr/bin/env python

"""
3.2-runDNAML.py

This script uses an iterative phylogenetic analysis to find sequences related
      to a set of known antibodies. Preprocessed sequences (including optional

Usage: 3.2-runDNAML.py -n native.fa -v germline_V
		       -i custom/input.phy
                       [-locus <H|K|L> -lib path/to/library.fa -f -h]


    Invoke with -h or --help to print this documentation.

    Required parameters:
    n		Fasta file containing the known sequences.
    v		Assigned germline V gene of known antibodes, for use in 
                   rooting the trees.

        * OR *

    i		Custom set of sequences to be analayzed, with natives and
                   germline sequence included, all ALIGNED and IN PHYLIP 
                   FORMAT. If -n and -v are supplied instead, sequences are
                   taken from output/sequences/ROOT-collected.fa and aligned 
                   automagically with MUSCLE. Manual alignment should be used
                   if ancestral sequences will be inferred.

    Optional parameters:	   
    locus	H (default): use V heavy / K: use V kappa / L: use V lambda
                   Ignored if the -lib option is used to supply a custom
                   library
    lib		Optional custom germline library (eg Rhesus or Mouse).

    Optional flags:
    f		Force a restart of the analysis, even if there are files from
                   a previous run in the working directory.

Created by Chaim A Schramm 2015-07-09.

Copyright (c) 2011-2015 Columbia University Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

from soanar import *
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline



def revertName(match):
    return lookup[ int(match.group(0)) - 1 ] #ids are 1-indexed, list is 0-indexed

def main():

    global inFile, lookup

    oldFiles = glob.glob("%s/infile"%prj_tree.phylo) + glob.glob("%s/outtree"%prj_tree.phylo) + glob.glob("%s/outfile"%prj_tree.phylo)
    if len(oldFiles) > 0:
        if force:
            for f in oldFiles:
                os.remove(f)
        else:
            sys.exit("Old files exist! Please use the -f flag to force overwrite.")
        

    if doAlign:

        #first create a working file to align and add the germline and natives
        shutil.copyfile("%s/%s-collected.fa"%(prj_tree.nt, prj_name), "%s/%s_to_align.fa"%(prj_tree.phylo, prj_name))
        handle = open( "%s/%s_to_align.fa"%(prj_tree.phylo, prj_name), "a" )
        handle.write( ">%s\n%s\n" % (germ_seq.id, germ_seq.seq) )
        for n in natives.values():
            handle.write( ">%s\n%s\n" % (n.id, n.seq) )
        handle.close()

        #now run muscle
        run_muscle            = MuscleCommandline( input="%s/%s_to_align.fa" % (prj_tree.phylo, prj_name), out="%s/%s_aligned.afa" % (prj_tree.phylo, prj_name) )
        run_muscle.maxiters   = 2
        run_muscle.diags      = True
        run_muscle.gapopen    = -5000.0 #code requires a float
        print run_muscle
        run_muscle()
        #thisVarHidesTheOutput = run_muscle()

        #change inFile variable so that remaining code is the same for both cases
        #It's probably really bad form to handle this in this way
        inFile = "%s/%s_aligned.afa" % (prj_tree.phylo, prj_name)


    #open the alignment to rename everything and find germline sequence
    #rename is to avoid possible errors with DNAML from sequence ids that are too long
    germ_pos = 1
    with open(inFile, "rU") as handle:
        if doAlign:
            aln = AlignIO.read(handle, "fasta")
        else: 
            try:
                aln = AlignIO.read(handle, "phylip")
            except:
                sys.exit("Please make sure custom input is aligned and in PHYLIP format")

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
    # O is outgroup root, followed by position of the germline in the alignment
    # 5 tells DNAML to do the ancestor inference
    # Y starts the run
    with open("%s/dnaml.in"%prj_tree.phylo, "w") as handle:
        seed = random.randint(0,1e10) * 2 + 1 #seed must be odd
        handle.write("J\n%d\n3\nO\n%d\n5\nY\n" % (seed, germ_pos))


    # change to work directory so DNAML finds "infile" and puts the output where we expect
    os.chdir(prj_tree.phylo)
    with open("%s/dnaml.in"%prj_tree.phylo, "rU") as pipe:
        subprocess.call([DNAML], stdin=pipe)

    #revert names in tree
    with open("%s/outtree"%prj_tree.phylo, "rU") as intree:
        mytree = intree.read()
    fixedtree = re.sub("\d{10}", revertName, mytree)
    with open("%s/%s.tree"%(prj_tree.out,prj_name), "w") as outtree:
        outtree.write(fixedtree)

    #revert names in out file
    with open("%s/outfile"%prj_tree.phylo, "rU") as instuff:
        mystuff = instuff.read()
    fixedstuff = re.sub("\d{10}", revertName, mystuff)
    with open("%s/%s.dnaml.out"%(prj_tree.logs,prj_name), "w") as outstuff:
        outstuff.write(fixedstuff)
        
    #clean up
    os.remove("infile")
    os.remove("outfile")
    os.remove("outtree")



if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)


	#check forcing parameter
	force = False
	flag = [x for x in ["f", "-f", "--f", "force", "-force", "--force"] if q(x)]
	if len(flag)>0:
		sys.argv.remove(flag[0])
		force = True


	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	= fullpath2last_folder(prj_tree.home)


	#get the parameters from the command line
	dict_args = processParas(sys.argv, n="natFile", v="germlineV", locus="locus", lib="library", i="inFile")
	defaults = dict(locus="H", library="")
	natFile, germlineV, locus, library, inFile = getParasWithDefaults(dict_args, defaults, "natFile", "germlineV", "locus", "library", "inFile")

        doAlign = True
	if natFile is None or germlineV is None:
            if inFile is None:
                print __doc__
		sys.exit(0)
            else:
                doAlign = False
        else:
            if inFile is not None:
                print __doc__
                sys.exit(0)
            else:
                #load native sequences
                natives = load_fastas(natFile)

                #load germline sequence
                if not os.path.isfile(library):
                    if locus in dict_vgerm_db.keys():
			library = dict_vgerm_db[locus]
                    else:
			print "Can't find custom V gene library file!"
			sys.exit(1)

                germ_dict = load_fastas(library)
                try:
                    germ_seq = germ_dict[germlineV]
                except:
                    print "Specified germline gene (%s) is not present in the %s library!\n" % (germlineV, library)
                    sys.exit(1)
		

	main()

