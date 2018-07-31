#!/usr/bin/env python

"""
5.1-repick_lineage_representative.py

This script uses the pseudo-lineages identified by 2.4-cluster_into_groups.py.
      It aligns each lineage using the full V(D)J sequence and selects a new
      lineage representative based on the sequences which is closest to the
      lineage consensus.

Usage: 5.1-repick_lineage_representative.py [ -i input.fa -o output.fa -m 2 ]

    All options are optional, see below for defaults.
    Invoke with -h or --help to print this documentation.

    i    File with lineage-anotated sequences to align. Default is output from 2.4:
             output/sequences/nucleotide/<project>_goodVJ_unique_lineageNotations.fa
    o    Where to save re-picked lineage representatives. Default is
             output/sequences/nucleotide/<project>_consensusLineageRepresentatives.fa
    m    Minimum number of sequences in psuedo-lineage to keep. Default is 2 (discard
             singletons).


Created by Chaim A Schramm on 2016-05-24.
Added to SONAR as part of mGSSP on 2017-02-24.
Modified to use VSearch instead of USearch by CAS on 2018-07-30.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from collections import OrderedDict
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline

try:
	from sonar.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/mGSSP")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *


def getHandle( item ):
    
    if item in handleList:
        return handleList[item]
    else:
        if len(handleList) >= 250: #being conservative here; increases runtime but decreases possibility of error
            i,h = handleList.popitem(False)
            h.close()
        handleList[item] = open("%s/%s.fa"%(prj_tree.lineage, item), "a")
        return handleList[item]


def main():

    #first, open the input file and parse into pseudo-lineages
    lineages = dict()
    count = 0
    for sequence in SeqIO.parse(open(inFile, "rU"), "fasta"):
        info = re.search(" size=(\d+) lineage_num=(\d+).*lineage_size=(\d+)", sequence.description)
        if info:
            if int(info.group(3)) >= minSeqs:
                if info.group(2) not in lineages:
                    lineages[info.group(2)] = dict( size=int(info.group(3)), desc=dict() )

                lineages[info.group(2)]['desc'][sequence.id] = sequence.description #because otherwise the alignment loses it

                #derep 'singletons' (seq size and lin size are the same) don't need to 
                #    be renamed for usearch and don't need saved handles
                if info.group(1) == info.group(3):
                    with open("%s/%s.fa"%(prj_tree.lineage ,info.group(2)), "w") as handle:
                        SeqIO.write([sequence], handle, 'fasta')
                else:
                    #add size notation in a way usearch can understand
                    sequence.id += ";size=%s" % info.group(1) #it's a string because that's what re.search returns
                    SeqIO.write([sequence], getHandle(info.group(2)), 'fasta')

        count += 1
        if count % 100000 == 0: print "Processed %d sequences in %d lineages so far..." % (count, len(lineages))

    #cleanup
    for h in handleList.itervalues(): h.close()
    handleList.clear()


    #go through each lineage and do the alignment
    reps = list()
    FNULL = open(os.devnull, 'w') #don't clutter up output with tons of usearch messages
    print "Starting usearch..."
    for lin in sorted(lineages, key=lambda num: lineages[num]['size'], reverse=True):

        #save time on singletons (if they weren't excluded by minSeq)
        if len(lineages[lin]['desc']) == 1:
            with open("%s/%s.fa" % (prj_tree.lineage, lin), "rU") as handle:
                reps.append( SeqIO.read(handle,'fasta') )
        else:
            #cluster and rapid align with vsearch
            subprocess.call([usearch, "-cluster_size", "%s/%s.fa" % (prj_tree.lineage, lin), 
                             "-id", "0.97", "-sizein", "-sizeout",
                             "-msaout", "%s/%s_msa.fa"%(prj_tree.lineage, lin), "-clusterout_sort"],
                             stdout=FNULL, stderr=subprocess.STDOUT)

            #extract biggest cluster
            with open("%s/%s_msa.fa"%(prj_tree.lineage, lin), "rU") as allClusters:
                with open("%s/%s_msaBiggest.fa"%(prj_tree.lineage, lin), "w") as biggestOnly:
                    blank = allClusters.next()
                    for line in allClusters:
                        if "consensus" in line:
                            break
                        biggestOnly.write(line)
            
            #open the msa
            with open("%s/%s_msaBiggest.fa" % (prj_tree.lineage, lin), "rU") as handle:
                aln = AlignIO.read(handle, "fasta")
                
            #add derep size to alignment as weighting
            for rec in aln:
                rec.annotations['weight'] = int( rec.id.split(";")[1].split("=")[1] )

            summary_align = AlignInfo.SummaryInfo(aln)
            pssm = summary_align.pos_specific_score_matrix()

            #score each sequence and save the best one
            scores = dict()
            for record in aln:
                myScore = 0
                for i,l in enumerate(record):
                    myScore += pssm[i][l]
                scores[record] = myScore
            d=sorted(aln, key=lambda rec: scores[rec], reverse=True) #reverse -> get max
            d[0].seq = d[0].seq.ungap("-") #remove gaps
            d[0].id = d[0].id.split(";")[0] #remove usearch size annotation
            d[0].id = re.sub("^\*","",d[0].id) #get rid of possible annotation from vsearch
            d[0].description = lineages[lin]['desc'][d[0].id] #restore original info
            reps.append( d[0] )
            
    
    #write output
    with open( outFile, "w" ) as handle:
        SeqIO.write( reps, handle, "fasta" )




if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)


	#save command line
	cmdLine = sys.argv

        prj_tree = ProjectFolders(os.getcwd())
        prj_name = fullpath2last_folder(prj_tree.home)


	# get parameters from input
	dict_args = processParas(sys.argv, i="inFile", o="outFile", m="minSeqs")
        default_args = dict(inFile="%s/%s_goodVJ_unique_lineageNotations.fa" % (prj_tree.nt, prj_name),
                           outFile="%s/%s_consensusLineageRepresentatives.fa" % (prj_tree.nt, prj_name),
                           minSeqs = 2)
	inFile, outFile, minSeqs = getParasWithDefaults(dict_args, default_args, "inFile", "outFile", "minSeqs")

        if not os.path.exists(inFile):
                print __doc__
		sys.exit(0)
        
        #since we're having to open/close/reopen filehandles to manage system limits
        #  clear out any old files so we're not double-counting data from multiple runs
        for f in glob.glob("%s/*fa"%prj_tree.lineage):
            os.remove(f)
        #make dir for msas but ignore error if it already exists
        try:
            os.mkdir("%s/msa"%prj_tree.lineage)
        except:
            pass

        handleList = OrderedDict()
        
	#log command line
	logCmdLine(sys.argv)	

	main()

