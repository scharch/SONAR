#!/usr/bin/env python3

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

Usage: 2.4-cluster_into_groups.py [ --id <90> --gaps <0> --natives FASTA (-v IGHV -j IGHJ) (--full FASTA --cdr3 FASTA) ]

Options:
    --id <90>         Clustering threshold (%) for CDR3 sequence identity (nucleotide).
                         [default: 90]
    --gaps <0>        Maximum number of (amino acid) in-dels to allow between CDR3
                         sequences in the same group. NOTE: This is implemented in
                         nucleotide space using vsearch's --maxgaps parameter, so does 
                         not guarantee an even codon in-del in the alignment used for 
                         clustering! [default: 0]
    --natives FASTA   Fasta file with nucleotide CDR3 sequences of known antibodies,
                         to be clustered together with the NGS data. Useful for focusing
                         on a known lineage.
    -v IGHV           V gene used by the known antibodies (without allele, eg: "IGHV1-2").
                         If not supplied, will be extracted from the fasta file by looking
                         for "v_call=IG[HKL]V..." in the fasta definition line. (The old
                         format of "V_gene=" will also still work.)
    -j IGHJ           J gene used by the known antibodies (without allele, eg: "IGHJ2").
                         If not supplied, will be extracted from the fasta file by looking
                         for "j_call=IG[HKL]J..." in the fasta definition line. (The old
                         format of "J_gene=" will also still work.)
    --full FASTA      Fasta file with full length sequences, if clustering something other
                         than SONAR's "goodVJ_unique" file. Must specify --cdr3, as well.
    --cdr3 FASTA      CDR3 nucleotide sequences extracted from custom sequence file if
                         not using SONAR's "goodVJ_unique" file. Must specify --full, too.


Created by Chaim A Schramm on 2015-04-27.
Edited by CAS 2017-07-26 to allow for automated extraction of V/J genes for native
                         antibodies and to exclude native-only clusters.
Edited by CAS 2018-06-13 to allow custom target file.
Modified to use VSearch by CAS 2018-07-30.
Changed regex for V genes and added special case handling for OR and (II)/(III) genes
                         by CAS 2018-08-28.
Edited to use Py3 and DocOpt by CAS 2018-08-29.
Updated for AIRR-format compatibility by CAS 2018-10-18.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
from collections import *
import airr

try:
	from sonar.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/lineage")
	sys.path.append(find_SONAR[0])
	from sonar.lineage import *


def main():

	#open logfile
	log = open("%s/cluster_into_groups.txt" % prj_tree.logs, "w")

	clusterLookup = dict()
	centroidData = dict()

	#first, open the input file and parse into groups with same V/J
	vj_partition = dict()
	cdr3_info = dict()
	seqSize = Counter()

	#start off by getting size annotations
	for read in generate_read_fasta(arguments['--full']):
		seqSize[read.id] = 1	
		check = re.search( "cluster_count=(\d+)", read.description)
		if check:
			seqSize[read.id] = int(check.group(1))

			
	gene_pat = re.compile("(?:v_call|V_gene)=IG([HKL]V[^*]+).*(?:j_call|J_gene)=IG([HKL]J\d)")
	for sequence in SeqIO.parse(open(arguments['--cdr3'], "rU"), "fasta"):
		genes = re.search(gene_pat, sequence.description)
		if genes:
			key = genes.group(1) + "_" + genes.group(2)
			key = re.sub("[()/]","",key) #so /OR or (II) genes don't screw up the file system
			if key not in vj_partition:
				temp = "%s/%s.fa"%(prj_tree.lineage, key) 
				vj_partition[key] = { 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

			vj_partition[key]['count'] += 1
			vj_partition[key]['ids'].append(sequence.id)
			cdr3_info[sequence.id] = { 'cdr3_len' : int(len(sequence.seq)/3), 'cdr3_seq' : sequence.seq.translate() }

			#make sizes available to vsearch
			sequence.id += ";size=%d" % seqSize[sequence.id] #do this even if there's no label
									 #so I don't need to divide the cases for vsearch
			#and write
			SeqIO.write([sequence], vj_partition[key]['handle'], 'fasta')
		else:
			print("Couldn't find V and J genes for %s %s, skipping..." % (sequence.id, sequence.description))


	natives = dict()
	if arguments['--natives'] is not None:
		natives = load_fastas(arguments['--natives'])
		for n, s in natives.items():
			if arguments['-v'] is not None:
				key = arguments['-v'] + "_" + arguments['-j']
			else:
				genes = re.search(gene_pat, s.description)
				if genes:
					key = genes.group(1) + "_" + genes.group(2)
				else:
					sys.exit("Can't find V and J gene annotations for native sequence %s. Please specify using the -v and -j parameters." % n)

			key = re.sub("[()/]","",key) #wouldn't expect this to be relevant for natives, but just in case...

			if key not in vj_partition:
				print( "No NGS sequences with the same V/J genes as native sequence %s (%s); skipping..." % (n, key) )
				continue

			seqSize[ n ] = 0
			s.id += ";size=1"
			vj_partition[key]['count'] += 1
			vj_partition[key]['ids'].append( n )
			cdr3_info[ n ] = { 'cdr3_len' : int(len(s.seq)/3), 'cdr3_seq' : s.seq.translate() }
			SeqIO.write([ s ], vj_partition[key]['handle'], 'fasta')

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
			if single in natives: centroidData[single]['nats'] = [single] #this shouldn't be possible, given we skipped natives with unique V/J combos above...
			clusterSizes[ single ] = seqSize[ single ]
			continue

		#cluster with vsearch
		subprocess.call([usearch, "-cluster_size", vj_partition[group]['file'], 
				 "-id", str(arguments['--id']/100.0), "-iddef", "1", #iddef=1 excludes gaps from also counting as mismatches
				 "-maxgaps", str(arguments['--gaps']*3),
				 "-sizein", "-uc", "%s/%s.uc"%(prj_tree.lineage, group),
				 "-leftjust", "-rightjust"], #left/right forces our pre-determined CDR3 borders to match 
				stdout=log, stderr=subprocess.STDOUT)

		#now reconstruct pseudo-lineages
		myGenes = group.split("_")
		with open("%s/%s.uc"%(prj_tree.lineage, group), "rU") as handle:
			uc = csv.reader( handle, delimiter=sep )
			for row in uc:
				#first get rid of size annotations
				hit  = re.sub(";size=\d+.*","",row[8])
				cent = re.sub(";size=\d+.*","",row[9]) # just a * for S rows, use hit as cent

				if cent in natives:
					#I am excluding clusters with natives as centroids because they only contain other natives
					#Since I've assigned them a size of 1, cluster_size will try them last, and since I added them
					# to the bottom of the file after all the NGs sequences, they should also be the last of any
					# other singletons to be tried, as well.
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
		writer.writerow([ "clone_id", "sequence_id", "v_call", "j_call", "junction_length_aa", 
				  "junction_aa", "clone_count", "included_mAbs" ])
		for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
			centroidData[centroid]['rank'] = rank+1
			writer.writerow([ "%05d"%(rank+1), centroid, centroidData[centroid]['vgene'], centroidData[centroid]['jgene'], 
					  cdr3_info[centroid]['cdr3_len'], cdr3_info[centroid]['cdr3_seq'], size, ",".join(centroidData[centroid]['nats']) ])

	#do sequence output
	notationFile = re.sub( "\.f.+", "_lineageNotations.fa", arguments['--full'] )
	repFile	     = re.sub( "\.f.+", "_lineageRepresentatives.fa", arguments['--full'] )
    
	rep_seqs = []
	with open( notationFile, "w" ) as handle:
		for read in generate_read_fasta(arguments['--full']):
			if ";" in read.id:
				read.id = read.id[0:8] #this is for raw USearch output with size annotations
						       #shouldn't be relevant in pipeline context
			if read.id not in clusterLookup: continue
			read.description += " clone_id=%05d clone_rep=%s clone_count=%d" % ( centroidData[clusterLookup[read.id]]['rank'], 
													   clusterLookup[read.id], clusterSizes[clusterLookup[read.id]] )
			SeqIO.write([read],handle,"fasta")
			if read.id in centroidData:
				rep_seqs.append(read)
					
	with open( repFile, "w" ) as handle:
		#use a sort to put them out in order of lineage rank (ie size)
		SeqIO.write( sorted(rep_seqs, key=lambda cent: centroidData[cent.id]['rank']), handle, "fasta" )

	#do AIRR output
	if os.path.dirname(arguments['--full']) == prj_tree.nt:
		if os.path.isfile("%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name)):
			withLin = airr.derive_rearrangement( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name),
							     fields=["clone_id", "clone_count"])
			for r in airr.read_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) ):
				if r['sequence_id'] in clusterLookup:
					r['clone_id']	 = "%05d"%centroidData[ clusterLookup[ r['sequence_id'] ] ][ 'rank' ]
					r['clone_count'] = clusterSizes[ clusterLookup[ r['sequence_id'] ] ]
				else:
					#prevent mix-and-match data if this gets run multiple times with multiple settings
					r['clone_id']	 = ""
					r['clone_count'] = ""
					
				withLin.write(r)
			withLin.close()
			os.rename( "updateRearrangements.tsv", "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) )

	log.close()



if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--id']  = int( arguments['--id'] )
	arguments['--gaps'] = int( arguments['--gaps'] )
	
	if arguments['--natives'] is not None:
		if os.path.isfile( arguments['--natives'] ):
			#working with known sequences, check for manual input of V and J genes
			if arguments['-v'] is None:
				print( "V and J genes of native sequences not specified; will try to parse from fasta file..." )
			elif not re.search("V", arguments['-v']) or not re.search("J", arguments['-j']):
				sys.exit("Cannot recognize input genes for known sequences.")
			elif re.search("\*", arguments['-v']) or re.search("\*", arguments['-j']):
				sys.exit("Please input germ line genes without allele (eg IGHV1-2 or IGKJ1)")
		else:
			sys.exit("Can't find native sequence file %s" % arguments['--natives'])

				
	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	if arguments['--full'] is not None:
		if not os.path.isfile(arguments['--full']):
			sys.exit("Cannot find sequence file %s" % arguments['--full'])
		if not os.path.isfile(arguments['--cdr3']):
			sys.exit("Cannot find CDR3 sequence file %s" % arguments['--cdr3'])
	else:
		arguments['--full'] = "%s/%s_goodVJ_unique.fa" % (prj_tree.nt, prj_name)
		arguments['--cdr3'] = "%s/%s_goodCDR3_unique.fa" % (prj_tree.nt, prj_name)

	#log command line
	logCmdLine(sys.argv)

	
	#create necessary directories if they don't exist
	os.makedirs(prj_tree.logs, exist_ok=True)
	os.makedirs(prj_tree.lineage, exist_ok=True)
	os.makedirs(prj_tree.tables, exist_ok=True)

    
	main()

