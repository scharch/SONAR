#!/usr/bin/env python

"""
2.4-cluster_into_groups.py

This script uses CDR3 identity to group unique sequences from a given data set
      into pseudo-lineages that can help define groups of related B cells.
      Uses output/sequences/nucleotide/<project>_goodCDR3_unique.fa as required
      input.
      The default threshold of 90% identity with no in-dels is probably useful
      for most cases, but more stringent or lenient criteria may be
      more appropriate in many cases.
      Sequences are first grouped by unique V and J gene assignments and then
      USearch is used to cluster the CDR3 sequences.

Usage: 2.4-cluster_into_groups.py [ -id 90 -gaps 0
                                    -n natives.fa -v nat_v_gene -j nat_j_gene
                                    -full sequences.fa -cdr3 cdr3_nt.fa ]

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
    natives.fa  Fasta file with nucleotide CDR3 sequences of known antibodies,
                   to be clustered together with the NGS data. Useful for focusing
                   on a known lineage.
    nat_v_gene  V gene used by the known antibodies (no allele, eg: "HV1-2").
                   If not supplied, will be extracted from the fasta file by
                   looking for "V_gene=IG[HKL]V..."
    nat_j_gene  J gene used by the known antibodies (no allele, eg: "HJ2")
                   If not supplied, will be extracted from the fasta file by
                   looking for "J_gene=IG[HKL]J..."
    full        Fasta file with full length sequences, if clustering something
                   other than unique ORFs. Must specify CDR3 file, as well.
    cdr3        CDR3 nucleotide sequences extracted from custom sequence file if
                   not using unique ORFs. Must specify -full as well.


Created by Chaim A Schramm on 2015-04-27.
Edited by CAS 2017-07-26 to allow for automated extraction of V/J genes for native
                         antibodies and to exclude native-only clusters.
Edited by CAS 2018-06-13 to allow custom target file.
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
    global natFile, nat_genes, gene_pat, idLevel, maxgaps, cdr3File, fullLengthFile

    #open logfile
    log = open("%s/cluster_into_groups.txt" % prj_tree.logs, "w")

    clusterLookup = dict()
    centroidData = dict()

    #first, open the input file and parse into groups with same V/J
    vj_partition = dict()
    cdr3_info = dict()
    seqSize = Counter()

    #start off by getting size annotations
    seqFile = "%s/%s_goodVJ_unique.fa" % (prj_tree.nt, prj_name)
    if fullLengthFile != "": seqFile = fullLengthFile
    for read in generate_read_fasta(seqFile):
        seqSize[read.id] = 1	
        check = re.search( " size=(\d+)", read.description)
        if check:
                seqSize[read.id] = int(check.group(1))
 



    inputFile = "%s/%s_goodCDR3_unique.fa" % (prj_tree.nt, prj_name)
    if cdr3File != "": inputFile = cdr3File
    
    for sequence in SeqIO.parse(open(inputFile, "rU"), "fasta"):
        genes = re.search(gene_pat, sequence.description)
        if genes:
            key = genes.group(1) + "_" + genes.group(2)
            if key not in vj_partition:
		temp = "%s/%s.fa"%(prj_tree.lineage, key) 
                vj_partition[key] = { 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

	    vj_partition[key]['count'] += 1
	    vj_partition[key]['ids'].append(sequence.id)
	    cdr3_info[sequence.id] = { 'cdr3_len' : len(sequence.seq)/3 - 2, 'cdr3_seq' : sequence.seq.translate() }

	    #make sizes available to usearch (using v9 formatting here)
	    sequence.description += ";size=%d;" % seqSize[sequence.description] #do this even if there's no label
	                                                      #so I don't need to divide the cases for usearch

	    #and write
	    SeqIO.write([sequence], vj_partition[key]['handle'], 'fasta')


    natives = dict()
    try:
        natives = load_fastas(natFile)
        for n, s in natives.items():
                key = nat_genes
                genes = re.search(gene_pat, s.description)
                if genes:
                        key = genes.group(1) + "_" + genes.group(2)

                if re.match("_", key) or re.search("_$",key):
                        print("V and J genes not found for native sequence %s, skipping it. Please check file or input manually."%n)
                        continue
                
                if key not in vj_partition:
                        print "No NGS sequences with the same V/J genes as native sequence %s (%s); skipping..." % (n, key)
                        continue

		seqSize[ n ] = 0
		s.description += ";size=0;"
		vj_partition[key]['count'] += 1
		vj_partition[key]['ids'].append( n )
		cdr3_info[ n ] = { 'cdr3_len' : len(s.seq)/3 - 2, 'cdr3_seq' : s.seq.translate() }
		SeqIO.write([ s ], vj_partition[key]['handle'], 'fasta')
    except IOError:
        pass


    #now go through and cluster each V/J grouping
    clusterSizes = Counter()
    for group in vj_partition:

	#close the file handle
        vj_partition[group]['handle'].close()

        #save a bit of time for obvious singletons
        if vj_partition[group]['count'] == 1:
            single = vj_partition[group]['ids'][0]
            myGenes = group.split("_")
            clusterLookup[ single ] = single
            centroidData[ single ] = dict( vgene = myGenes[0], jgene = myGenes[1], nats=[] )
	    if single in natives: centroidData[single]['nats'] = [single]
	    clusterSizes[ single ] = seqSize[ single ]
            continue

        #cluster with usearch
        subprocess.call([usearch, "-cluster_fast", vj_partition[group]['file'], 
                         "-id", str(idLevel/100.0), "-maxgaps", str(maxgaps*3), "-sizein",
                         "-sort", "size", "-uc", "%s/%s.uc"%(prj_tree.lineage, group),
                         "-leftjust", "-rightjust"], #left/right forces our pre-determined CDR3 borders to match 
                        stdout=log, stderr=subprocess.STDOUT)

        #now reconstruct pseudo-lineages
        myGenes = group.split("_")
	with open("%s/%s.uc"%(prj_tree.lineage, group), "rU") as handle:
		uc = csv.reader( handle, delimiter=sep )
		for row in uc:
			#first get rid of size annotations and fasta def line (which is included as of usearch9)
			#hit  = re.sub(";size=\d+;.*","",row[8])
			#cent = re.sub(";size=\d+;.*","",row[9]) # just a * for S rows, use hit as cent
			hit  = re.sub(" .*","",row[8])
			cent = re.sub(" .*","",row[9]) # just a * for S rows, use hit as cent

                        if cent in natives:
                                continue
			elif row[0] == "S":
				if hit in natives:
                                        continue
				centroidData[ hit ] = dict( vgene = myGenes[0], jgene = myGenes[1], nats=[] )
			        clusterLookup[ hit ] = hit
				clusterSizes[ hit ] = seqSize[ hit ]
			elif row[0] == "H":
				clusterLookup[ hit ] = cent
				clusterSizes[ cent ] += seqSize[ hit ]
				if hit in natives:
					centroidData[ cent ][ 'nats' ].append( hit )
			else:
				break #skip "C" lines
        

    #now process all clusters and do tabular output
    with open( "%s/%s_lineages.txt" % (prj_tree.tables, prj_name), "w" ) as handle:
        writer = csv.writer(handle, delimiter=sep)
        writer.writerow([ "lineage_ID", "rep_seq_ID", "V_gene", "J_gene", "cdr3_len", 
                       "cdr3_aa_seq", "size", "included_mAbs" ])
        for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
            centroidData[centroid]['rank'] = rank+1
            writer.writerow([ "%05d"%(rank+1), centroid, centroidData[centroid]['vgene'], centroidData[centroid]['jgene'], 
                              cdr3_info[centroid]['cdr3_len'], cdr3_info[centroid]['cdr3_seq'], size, ",".join(centroidData[centroid]['nats']) ])

    #do sequence output
    notationFile = re.sub( "\.f.+", "_lineageNotations.fa", seqFile )
    repFile      = re.sub( "\.f.+", "_lineageRepresentatives.fa", seqFile )
    
    rep_seqs = []
    with open( notationFile, "w" ) as handle:
        for read in generate_read_fasta(seqFile):
            if ";" in read.id:
                read.id = read.id[0:8] #this is for raw USearch output with size annotations
                                       #shouldn't be relevant in pipeline context
	    if read.id not in clusterLookup: continue
            read.description += " lineage_num=%05d lineage_rep=%s lineage_size=%d" % ( centroidData[clusterLookup[read.id]]['rank'], 
                                                                                       clusterLookup[read.id], clusterSizes[clusterLookup[read.id]] )
            SeqIO.write([read],handle,"fasta")
            if read.id in centroidData:
                rep_seqs.append(read)

    with open( repFile, "w" ) as handle:
        #use a sort to put them out in order of lineage rank (ie size)
        SeqIO.write( sorted(rep_seqs, key=lambda cent: centroidData[cent.id]['rank']), handle, "fasta" )


    log.close()

if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
		print __doc__
		sys.exit(0)

	#log command line
	logCmdLine(sys.argv)

	# get parameters from input
	dict_args = processParas(sys.argv, id="idLevel", gaps="maxgaps", n="natFile", v="natV", j="natJ", cdr3="cdr3File", full="fullLengthFile")
	idLevel, maxgaps, natFile, natV, natJ, cdr3File, fullLengthFile = getParasWithDefaults(dict_args, dict(idLevel=90, maxgaps=0, natFile="", natV="", natJ="", cdr3File="", fullLengthFile=""),
                                                                     "idLevel", "maxgaps", "natFile", "natV", "natJ", "cdr3File", "fullLengthFile")

        if os.path.isfile(natFile):
            #working with known sequences, check for manual input of V and J genes
            if natV=="" or natJ=="":
                print "V and/or J gene of native sequences not specified; will try to parse from fasta file..."
            elif not re.search("V", natV) or not re.search("J", natJ):
                sys.exit("Cannot recognize input genes for known sequences.")
            elif re.search("\*", natV) or re.search("\*", natJ):
                sys.exit("Please input germ line genes without allele (eg HV1-2 or KJ1)")


        if (cdr3File=="" and fullLengthFile!="") or (cdr3File!="" and fullLengthFile==""):
                sys.exit("Please specify files for both full length sequences and extracted CDR3s.")
                

        nat_genes = natV + "_" + natJ
        #gene_pat = re.compile("V_gene=IG([HKL]V\d+(?:D|\/OR\d*)?Y?-(?:NL)?\d+).*J_gene=IG([HKL]J\d)")
        gene_pat = re.compile("V_gene=IG([HKL]V[^*]+).*J_gene=IG([HKL]J\d)")


        prj_tree = ProjectFolders(os.getcwd())
        prj_name = fullpath2last_folder(prj_tree.home)


	main()

