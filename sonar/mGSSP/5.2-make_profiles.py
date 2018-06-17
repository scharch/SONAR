#!/usr/bin/env python

"""

5.2-make_profiles.py

This script takes input sequences and generates germline mutability profiles
      for downstream analysis. The recommended workflow is to run the annotation
      module (1.1 through 1.4) followed by 2.4-cluster_into_groups.py and (optionally)
      5.1-repick_lineage_representatives.py. Then check sequences for frameshifts
      using utilities/checkForFrameshift.py, translate the sequences (utilties/quickTranslate.py)
      and run 2.1-calculate_id-div.pl on the translated sequences to identify and
      remove those with no amino acid substitutions. Then run this program.

Usage: 5.2-make_profiles.py <sequences.fa> [ -o profiles.txt -n 50 -g germV.fa -a -u ]

Options:
   -h --help                    Show this documentation
   <sequences.fa>               Processed ngs sequences to be used to build the profiles
   -o --output profiles.txt     Where to save output [default: profiles.txt]
   -n --numPerProfile 50        Number of reads to use in building each profile [default: 50]
   -g --germline germV.fa       Location of germline V sequences to use for building
                                   profiles. Expected as trimmed/padded AA sequences.
                                   [default: /ifs/scratch/c2b2/ls_nih/cs3037/spectrum/nonredundant_AA_alleles_includingNovel_20160609.fa]
   -a                           Input sequences are amino acid (don't translate) [default: False]
   -u --use-all                 Don't subset or run multiple profiles. Still uses the numPerProfile
                                   parameter as a minimum [default: False]

Created by Chaim Schramm on 2016-05-27.
Added to SONAR as part of mGSSP on 2017-02-24.
Copyright (c) 2011-2017 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys, numpy, re, csv, glob, os
from collections import defaultdict, Counter
from docopt import docopt
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

try:
	from sonar.mGSSP import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/mGSSP")
	sys.path.append(find_SONAR[0])
	from sonar.annotate import *

        
def main():
    
    #load sequences
    masterList = defaultdict( list )
    with open(arguments["<sequences.fa>"], 'rU') as handle:
        for sequence in SeqIO.parse(handle, "fasta"):

            #start with a special case where IMGT allele is misnamed
            sequence.description = re.sub("IGHV4-4\*0[78]", "IGHV4-59*11", sequence.description)

            #collapse distal alleles
            sequence.description = re.sub("(V\d-\d+)D", r'\1', sequence.description)

            gene = re.search("(IG[HKL]V\d-[^*]+)", sequence.description) #breaks all the ORs - don't care
            if gene:
                if not arguments["-a"]:
                    sequence.seq = sequence.seq.translate() #don't care about anything else in this script
                masterList[ gene.group(1) ].append( sequence )

    #replace list with array
    #have to build manually because otherwise numpy is turning SeqRecords
    #            into lists of chars (AAs), which causes random.choice
    #            to throw an error (non 1-D array)
    # (weird footnote: this _doesn't_ happen if 1 or more sequences in the
    #            don't have an even number of codons)
    # Anyway, only fix I can come up with is to manually place each SeqRecord
    #            in the array. We have to do it here, afterword, because until
    #            we've finished loading the sequences, I don't know how many
    #            there will be of each germline and numpy arrays have to be
    #            pre-allocated.

    for v in masterList.keys():
        a = numpy.empty( len(masterList[v]), dtype=object )
        for i in range(len(masterList[v])):
            a[i] = masterList[v][i]
        masterList[v] = a


    #load germlines
    germList = defaultdict( list )
    with open(arguments["--germline"], 'rU') as handle:
        for sequence in SeqIO.parse(handle, "fasta"):
            #collapse distal alleles and remove allele designation
            gene = re.sub("(V\d-\d+)D?\*.*", r'\1', sequence.id)
            germList[ gene ].append( sequence )
    
    
    #start output file
    outHandle = open(arguments["--output"], "w")
    output = csv.writer(outHandle, delimiter="\t")
    output.writerow( ["Vgene", "prof#", "pos", "germ", "freq"] + aa_list )

    #now let's start building profiles
    for v in sorted(masterList.keys()):

        if len(masterList[v]) < arguments["--numPerProfile"]:
            print "Skipping %s, not enough sequences..." % v
            continue #not enough data for a profile

        if v not in germList:
            print "Skipping %s, it's not in the germline database..." %v
            continue


        #we're going to take random overlapping subsets; how many is reasonable?
        #   the way I am doing this is completely arbitrary
        #   every 10 sequences, we'll add another iteration (or =20% for larger n)
        #   But cap it at 5, no need to get carried away
        numProfiles = int(numpy.ceil( (len(masterList[v])-arguments["--numPerProfile"])/(float(arguments["--numPerProfile"])/5) ))
        if numProfiles > 5: numProfiles = 5
        if arguments["--use-all"]:
            numProfiles = 1

        success = 0
        for i in range(numProfiles):


            seqs = [] + germList[v] #force a copy rather than an alias
            if arguments["--use-all"]:
                seqs += list(masterList[v])
                #with open("../log.txt", "a") as handle:
                #    handle.write( "\t".join([ arguments["<sequences.fa>"], v, str(len(masterList[v])) ]) + "\n" )
            else:
                #get our sequence subset, add the germlines, and write them
                #   to a temporary file for alignment
                seqs += list(numpy.random.choice(masterList[v], size=arguments["--numPerProfile"], replace=False))

            tempFile = "".join(numpy.random.choice(list("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"), size=8))
            with open("%s.fa"%tempFile, "w") as temp:
                SeqIO.write(seqs,temp,"fasta")

            clustal_cline = MuscleCommandline(input="%s.fa"%tempFile, out="%s.aln"%tempFile)
            try:
                stdout, stderr = clustal_cline()
            except:
                print "Error in alignment #%d for %s (skipping)" % (i+1, v)
                for f in glob.glob("%s.*"%tempFile): 
                    os.remove(f)
                continue

            alignment = AlignIO.read("%s.aln"%tempFile, "fasta")#"clustal")
            success += 1

            #Input order is not maintained, so we need a little
            #   kludge to find a germline sequences. Use the 
            #   first one to remove any insertions from the alignment
            germRow = 0
            for n, rec in enumerate(alignment):
                if rec.id in [g.id for g in germList[v]]:
                    germRow = n
                    break

            #look for gaps one at a time so we don't get tripped up by shifting indices
            gap = re.search( "-+", str(alignment[germRow].seq) )
            while (gap):
                alignment = alignment[:, 0:gap.start()] + alignment[:, gap.end():]
                gap = re.search( "-+", str(alignment[germRow].seq) )
                
            #Now we get BioPython to make a PSSM for us. To convert that into
            #    a mutability profile, we will delete the germline residue[s]
            #    at each position (but save what they were)
            germRes = defaultdict(Counter)
            summary_align = AlignInfo.SummaryInfo(alignment)
            pssm = summary_align.pos_specific_score_matrix(chars_to_ignore=['-','X'])
            for germ in germList[v]:
                for pos, residue in enumerate(germ):
                    if residue == "X":
                        continue
                    germRes[pos][residue] += 1
                    pssm[pos][residue] = 0

            #normalize and save
            numInProfile = arguments["--numPerProfile"]
            if arguments["--use-all"]:
                numInProfile = len(masterList[v])
            for p, pos in enumerate(pssm):
                germAA = ",".join([ x[0] for x in germRes[p].most_common() ])
                output.writerow( [ v, i+1, p+1, germAA, "%.5f"%(sum(pos.values())/numInProfile) ] + [ "%.5f"%(pos.get(r,0)/sum(pos.values())) if sum(pos.values()) > 0 else "0.00" for r in aa_list ] )
            
            #clean up
            for f in glob.glob("%s.*"%tempFile): 
                os.remove(f)

        print "Successfully built %d/%d profiles for %s!" % (success, numProfiles, v)


    outHandle.close()



if __name__ == '__main__':

    aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    #save command line
    cmdLine = sys.argv

    arguments = docopt(__doc__)
    arguments['--numPerProfile'] = int(arguments['--numPerProfile'])

    #log command line
    logCmdLine(sys.argv)
    
    main()

