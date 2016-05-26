#!/usr/bin/env python

"""
2.4-cluster_into_groups.py

This script uses CDR3 identity to group unique sequences from a given data set
      into pseudo-lineages that can help define groups of related B cells.
      The default threshold of 90% identity with no in-dels is probably useful
      for most cases, but more stringent or lenient criteria may be
      more appropriate in many cases.
      Sequences are first grouped by unique V and J gene assignments and then
      USearch is used to cluster the CDR3 sequences.

Usage: 2.4-cluster_into_groups.py [ -id 90 -gaps 0
                                    -n natives.fa -v nat_v_gene -j nat_j_gene ]

    All options are optional, see below for defaults.
    Invoke with -h or --help to print this documentation.

    id          Clustering threshold (%) for CDR3 sequence identity (nucleotide).
                   Default = 90.
    gaps        Maximum number of (amino acid) in-dels to allow between CDR3
                   sequences in the same group. Default = 0. NOTE: This is
                   implemented in nucleotide space using the USearch "maxgaps"
                   parameter, so does not guarantee an even codon in-del in the
                   alignment used for clustering! Also, usearch will still count
                   gaps as mismatches, so a CDR3 of 20AA will be counted as 95%
                   id to an identical-other-than-deletion 19AA CDR3. Set your
                   thresholds accordingly.
    natives.fa  Fasta file of known antibody sequences to be clustered together
                   with the NGS data. Useful for focusing on a known lineage.
    nat_v_gene  V gene used by the known antibodies (no allele, eg: "HV1-2")
    nat_j_gene  J gene used by the known antibodies (no allele, eg: "HJ2")


Created by Chaim A Schramm on 2015-04-27.
Copyright (c) 2011-2016 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from collections import *

try:
	from sonar.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/lineage")
	sys.path.append(find_SONAR[0])
	from sonar.lineage import *


def main():

    #need both CDR3 and full length - infile maybe not parameter
    global natFile, nat_genes, gene_pat, idLevel, maxgaps

    #open logfile
    log = open("%s/cluster_into_groups.txt" % prj_tree.logs, "w")

    clusterLookup = dict()
    centroidData = dict()

    #first, open the input file and parse into groups with same V/J
    vj_partition = dict()
    seqSize = Counter()
    for sequence in SeqIO.parse(open("%s/%s_goodCDR3_unique.fa" % (prj_tree.nt, prj_name), "rU"), "fasta"):
        genes = re.search(gene_pat, sequence.description)
        if genes:
            key = genes.group(1) + "_" + genes.group(2)
            if key in vj_partition:
                vj_partition[key][sequence.id] = sequence
            else:
                vj_partition[key] = { sequence.id : sequence }
	    #add sizes
	    seqSize[sequence.id] = 1	
	    check = re.search( " size=(\d+)", sequence.description)
	    if check:
		    seqSize[sequence.id] = check.group(1)


    natives = dict()
    try:
        natives = load_fastas(natFile)
        for n, s in natives.items():
            try:
                vj_partition[nat_genes][n] = s
            except KeyError:
                vj_partition[nat_genes] = { n : s }
	    seqSize[ n ] = 1
    except IOError:
        pass


    #now go through and cluster each V/J grouping
    clusterSizes = Counter()
    for group in vj_partition:

        #save a bit of time for obvious singletons
        if len(vj_partition[group]) == 1:
            single = vj_partition[group].values()[0]
            myGenes = group.split("_")
            clusterLookup[ single.id ] = single.id
            centroidData[ single.id ] = dict( vgene = myGenes[0], jgene = myGenes[1], cdr3_len = len(single.seq)/3 - 2, 
                                                cdr3_seq = single.seq.translate(), nats=[] )
	    clusterSizes[ single.id ] = seqSize[ single.id ]
            continue

        tempFile = open("%s/%s.fa"%(prj_tree.lineage, group), "w")
        SeqIO.write(vj_partition[group].values(), tempFile, "fasta")
        tempFile.close()

        #cluster with usearch
        subprocess.call([usearch, "-cluster_fast", "%s/%s.fa"%(prj_tree.lineage, group), 
                         "-id", str(idLevel/100.0), "-maxgaps", str(maxgaps),
                         "-sort", "size", "-uc", "%s/%s.uc"%(prj_tree.lineage, group),
                         "-leftjust", "-rightjust"], #left/right forces our pre-determined CDR3 borders to match 
                        stdout=log, stderr=subprocess.STDOUT)

        #now reconstruct pseudo-lineages
        myGenes = group.split("_")
        uc = csv.reader( open("%s/%s.uc"%(prj_tree.lineage, group), "rU"), delimiter=sep )
        for row in uc:
            if row[0] == "S":
                seq = vj_partition[group][row[8]]
                centroidData[ row[8] ] = dict( vgene = myGenes[0], jgene = myGenes[1], cdr3_len = len(seq.seq)/3 - 2, 
                                                cdr3_seq = seq.seq.translate(), nats=[] )
                clusterLookup[ row[8] ] = row[8]
		clusterSizes[ row[8] ] = seqSize[ row[8] ]
                if row[8] in natives:
                    centroidData[ row[8] ][ 'nats' ].append( row[8] )
            elif row[0] == "H":
                clusterLookup[ row[8] ] = row[9]
		clusterSizes[ row[9] ] += seqSize[ row[8] ]
                if row[8] in natives:
                    centroidData[ row[9] ][ 'nats' ].append( row[8] )
            else:
                break #skip "C" lines
        

    #now process all clusters and do tabular output
    with open( "%s/%s_lineages.txt" % (prj_tree.tables, prj_name), "w" ) as handle:
        writer = csv.writer(handle, delimiter=sep)
        writer.writerow([ "lineage_ID", "rep_seq_ID", "V_gene", "J_gene", "cdr3_len", 
                       "cdr3_aa_seq", "size", "included_mAbs" ])
        for rank, [centroid, size] in enumerate(clusterSizes.most_common()):
            centroidData[centroid]['rank'] = rank+1
            centroidData[centroid]['size'] = size
            tempDict = centroidData[centroid]
            writer.writerow([ rank+1, centroid, tempDict['vgene'], tempDict['jgene'], 
                           tempDict['cdr3_len'], tempDict['cdr3_seq'], size, ",".join(tempDict['nats']) ])

    #do sequence output
    rep_seqs = []
    with open( "%s/%s_goodVJ_unique_lineageNotations.fa" % (prj_tree.nt, prj_name), "w" ) as handle:
        for read in generate_read_fasta("%s/%s_goodVJ_unique.fa" % (prj_tree.nt, prj_name)):
            if ";" in read.id:
                read.id = read.id[0,8] #this is for raw USearch output with size annotations
                                       #shouldn't be relevant in pipeline context
	    if read.id not in clusterLookup: continue
            read.description += " lineage_num=%05d lineage_rep=%s lineage_size=%d" % ( centroidData[clusterLookup[read.id]]['rank'], 
                                                                                       clusterLookup[read.id], centroidData[clusterLookup[read.id]]['size'] )
            SeqIO.write([read],handle,"fasta")
            if read.id in centroidData:
                rep_seqs.append(read)

    with open( "%s/%s_lineageRepresentatives.fa" % (prj_tree.nt, prj_name), "w" ) as handle:
        #use a sort to put them out in order of lineage rank (ie size)
        SeqIO.write( sorted(rep_seqs, key=lambda cent: centroidData[cent.id]['rank']), handle, "fasta" )


    log.close()

if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, id="idLevel", gaps="maxgaps", n="natFile", v="natV", j="natJ")
	idLevel, maxgaps, natFile, natV, natJ = getParasWithDefaults(dict_args, dict(idLevel=90, maxgaps=0, natFile="", natV="", natJ=""),
                                                                     "idLevel", "maxgaps", "natFile", "natV", "natJ")

        if os.path.isfile(natFile):
            #ok working with known sequences, make sure V and J were input properly
            if not re.search("V", natV) or not re.search("J", natJ):
                sys.exit("Cannot recognize input genes for known sequences.")
            elif re.search("\*", natV) or re.search("\*", natJ):
                sys.exit("Please input germ line genes without allele (eg HV1-2 or KJ1)")
        

        nat_genes = natV + "_" + natJ
        gene_pat = re.compile("([HKL]V\d-[^*]+).*([HKL]J\d)")


        prj_tree = ProjectFolders(os.getcwd())
        prj_name = fullpath2last_folder(prj_tree.home)


	main()

