#!/usr/bin/env python3

"""
2.4-cluster_into_groups.py

This script uses CDR3 identity to group unique sequences from a given data set
      into pseudo-lineages that can help define groups of related B cells.
      Uses output/sequences/nucleotide/<project>_goodCDR3_unique.fa and
      output/sequences/nucleotide/<project>_goodVJ_unique.fa as default input.

      The default threshold of 90% identity with no in-dels is probably useful
      for most cases, but more stringent or lenient criteria may sometimes be
      more appropriate.

      Sequences are first grouped by unique V and J gene assignments and then
      VSearch is used to cluster the CDR3 sequences.

Usage: 2.4-cluster_into_groups.py [ --id <90> --gaps <0> --natives FASTA (-v IGHV -j IGHJ) (--full FASTA --cdr3 FASTA) --singlecell -t 1 ]

Options:
    --id <90>         Clustering threshold (%) for CDR3 sequence identity (nucleotide).
                         [default: 90]
    --gaps <0>        Maximum number of in-dels to allow between CDR3 sequences in the
                         same group. NOTE: This is implemented in nucleotide space using
                         vsearch's --maxgaps parameter, which counts gap openings, rather
                         than gap columns, so probably shouldn't be set to more than 1.
                         There is no guarantee an even codon in-del in the alignment used for
                         for clustering! Also, vsearch will still count the gap columns as
                         mismatches, so a CDR3 of 20AA will be counted as 95% id to an
                         identical-other-than-deletion 19AA CDR3. Set your --id
                         threshold accordingly. [default: 0]
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
    --singlecell      A flag to indicate single cell data - heavy and light chain data will
                         be used jointly to define clones. output/tables/rearrangements_single-cell.tsv
                         (output of 1.5) must be present. [default: False]
    -t 1              Number of threads used [default: 1]


Created by Chaim A Schramm on 2015-04-27.
Edited by CAS 2017-07-26 to allow for automated extraction of V/J genes for native
                         antibodies and to exclude native-only clusters.
Edited by CAS 2018-06-13 to allow custom target file.
Modified to use VSearch by CAS 2018-07-30.
Changed regex for V genes and added special case handling for OR and (II)/(III) genes
                         by CAS 2018-08-28.
Edited to use Py3 and DocOpt by CAS 2018-08-29.
Updated for AIRR-format compatibility by CAS 2018-10-18.
Added multithreading by CAS 2019-02-05.
Generalized germline gene regexes by CAS 2020-01-02.
Added joint heavy-light clonal partitioning for single cells by CAS 2020-01-16.

Copyright (c) 2011-2020 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
from collections import *
from multiprocessing import Pool
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import airr

try:
	from SONAR.lineage import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/lineage")
	sys.path.append(find_SONAR[0])
	from SONAR.lineage import *


# a utility function to get us a slice of an iterator, as an iterator
# when working with iterators maximum lazyness is preferred
# from https://stackoverflow.com/a/44502827
def iterator_slice(iterator, length):
	iterator = iter(iterator)
	slice = 0
	while True:
		res = tuple(itertools.islice(iterator, length))
		if not res:
			break
		slice += 1
		yield slice, res



def processClusters( iter_tuple ):

	count, chunk = iter_tuple
	print("Processing chunk #%d..."%count)

	global natives

	clusterLookup = dict()
	centroidData = dict()
	clusterSizes = Counter()

	for cluster in chunk:

		#save a bit of time for obvious singletons
		if cluster['count'] == 1:
			single = cluster['ids'][0]
			myGenes = cluster['group'].split("_")
			clusterLookup[ single ] = single
			centroidData[ single ] = dict( vgene = myGenes[0], jgene = myGenes[1], nats=[] )
			if single in natives: centroidData[single]['nats'] = [single] #this shouldn't be possible, given we skipped natives with unique V/J combos above...
			clusterSizes[ single ] = 1#seqSize[ single ]
			continue

		#cluster with vsearch
		subprocess.call([vsearch, "-cluster_size", cluster['file'],
				 "-id", str(arguments['--id']/100.0),
				 "-maxgaps", str(arguments['--gaps']),
				 "-sizein", "-uc", "%s/%s.uc"%(prj_tree.lineage, cluster['group']),
				 "-minseqlength", "15", #lets us capture CDR3s down to 3 aa
				 "-leftjust", "-rightjust", #left/right forces our pre-determined CDR3 borders to match
				 "-quiet"] #supress screen clutter
				)

		#now reconstruct pseudo-lineages
		myGenes = cluster['group'].split("_")
		with open("%s/%s.uc"%(prj_tree.lineage, cluster['group']), "r") as handle:
			uc = csv.reader( handle, delimiter=sep )
			for row in uc:
				#first get rid of size annotations
				hit  = re.sub(";size=\d+.*","",row[8])
				cent = re.sub(";size=\d+.*","",row[9]) # just a * for S rows, use hit as cent

				if row[0] == "S":
					centroidData[ hit ] = dict( vgene = myGenes[0], jgene = myGenes[1], nats=[] )
					if hit in natives:
						centroidData[ hit ][ 'nats' ].append(hit)
					else:
						clusterSizes[ hit ] = 1#seqSize[ hit ]
					clusterLookup[ hit ] = hit
				elif row[0] == "H":
					clusterLookup[ hit ] = cent
					if hit in natives:
						centroidData[ cent ][ 'nats' ].append( hit )
					else:
						clusterSizes[ cent ] += 1#seqSize[ hit ]
				else:
					break #skip "C" lines

	return { 'cl':clusterLookup, 'cd':centroidData, 'cs':clusterSizes }



def jointClonality(clusters, cells, cdr3Info, natCells):

	#start by creating a graph connecting each pair of 'chain clones' that have been seen
	# together in an individual cell
	# use the `cluster` dictionary to find the centroid to which each chain has been assigned
	cloneGraph = Graph()
	for c in cells:
		for l1 in range(len(cells[c])):
			for l2 in range(l1+1, len(cells[c])):
				cloneGraph.add_edge( clusters[ cells[c][l1] ], clusters[ cells[c][l2] ] )


	#now decompose the graph into "maximal cliques" to identify 'cell clones'
	# probably overkill since the graph will mostly be disjoint, with clearly identifiable
	# clones, but this way I don't have to think about edge cases (no pun intended)
	cellClones = Counter()
	for cliq in find_cliques(cloneGraph):
		cellClones[ frozenset(cliq) ] = -1


	#Go back through the list of cells and assign them to clones
	assignments = dict()
	natClones = defaultdict(list)
	ambiguous = 0
	for c in cells:
		possibleClones = []
		for cc in cellClones:
			if frozenset([clusters[s] for s in cells[c]]) <= cc:
				possibleClones.append(cc)
		if len(possibleClones) == 1:
				assignments[c] = possibleClones[0]
				if cellClones[ possibleClones[0] ] < 0:
					cellClones[ possibleClones[0] ] = 0 #kludge to distinguish between those unseen and those with only "natives"
				if c in natCells:
					natClones[ possibleClones[0] ].append( c )
				else:
					cellClones[ possibleClones[0] ] += 1
		else:
			#no matches found, probably ambiguous
			# but check if its a singleton that never appeared with any other chain
			if len(possibleClones) == 0 and len(cells[c]) == 1:
				assignments[c] = frozenset([clusters[s] for s in cells[c]])
				if c in natCells:
					natClones[ frozenset([clusters[s] for s in cells[c]]) ].append( c )
					cellClones[ frozenset([clusters[s] for s in cells[c]]) ] = 0
				else:
					cellClones[ frozenset([clusters[s] for s in cells[c]]) ] = 1

			else:
				#no assignment
				ambiguous += 1


	#issue warning if too many are ambiguous/unassignable
	if ambiguous > len(cells)/20:
		print( "Warning: More than 5% of cells had ambiguous or unassignable clonality.", file=sys.stderr)

	#finally, collect the detailed 'chain clone' data for each 'cell clone'
	cloneInfo = dict()
	for cc in cellClones:
		nt = []
		aa = []
		for chain in cc:
			nt.append( f"{cdr3Info[chain]['genes']}:{cdr3Info[chain]['cdr3_seq']}" )
			aa.append( f"{cdr3Info[chain]['genes']}:{cdr3Info[chain]['cdr3_seq'].translate()}" )

		cloneInfo[cc] = { 'aa':",".join(sorted([x.upper() for x in aa])), 'nt':",".join(sorted([x.upper() for x in nt])), 'nats':"" }
		if cc in natClones:
			cloneInfo[cc]['nats'] = natClones[cc]

	return assignments, cloneInfo, cellClones



def main():

	#first, open the input file and parse into groups with same V/J
	vj_partition = dict()
	cdr3_info = dict()
	seqSize = Counter()

	#AIRR TSV has to be parsed differently
	if arguments['--singlecell']:

		cell_dict = defaultdict(list)
		for r in filterAirrTsv(f"{prj_tree.tables}/{prj_name}_rearrangements_single-cell.tsv", [{'column':'cell_status','list':["^(?:(?!multi|none).)*$"]}]):

			#short CDR3s can't be clustered and probably indicate a bad sequence anyway
			if r['junction_length'] < 15:
				continue

			#skip sequences with missing gene assignments
			if r['v_call'] == "" or r['j_call'] == "":
				continue

			#for network analysis
			cell_dict[ r['cell_id'] ].append(r['sequence_id'])

			#get size
			seqSize[ r['sequence_id'] ] = 1
			if not r['duplicate_count'] == "":
				seqSize[ r['sequence_id'] ] = r['duplicate_count']

			#get gene assignments
			key = r['v_call'].split("*")[0] + "_" + r['j_call'].split("*")[0]
			if key not in vj_partition:
				temp = "%s/%s.fa"%(prj_tree.lineage, key)
				vj_partition[key] = { 'group':key, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

			vj_partition[key]['count'] += 1
			vj_partition[key]['ids'].append(r['sequence_id'])
			cdr3_info[ r['sequence_id'] ] = { 'genes' : key, 'cdr3_seq' : Seq.Seq(r['junction']) }

			#create a sequence object
			tempSeq = SeqRecord( id=r['sequence_id'], seq=Seq.Seq(re.sub("[-.+]","",r['junction'])) )
			tempSeq.id += ";size=%d" % seqSize[ r['sequence_id'] ] #do this even if there's no label
									                               #so I don't need to divide the cases for vsearch
			#and write
			SeqIO.write([tempSeq], vj_partition[key]['handle'], 'fasta')

	else:
		#start off by getting size annotations
		for read in generate_read_fasta(arguments['--full']):
			seqSize[read.id] = 1
			check = re.search( "cluster_count=(\d+)", read.description)
			if check:
				seqSize[read.id] = int(check.group(1))

		v_gene_pat = re.compile("(?:v_call|V_gene)=([^*\s]+)")
		j_gene_pat = re.compile("(?:j_call|J_gene)=([^*\s]+)")
		for sequence in SeqIO.parse(open(arguments['--cdr3'], "r"), "fasta"):
			vgene = re.search(v_gene_pat, sequence.description)
			jgene = re.search(j_gene_pat, sequence.description)
			if vgene and jgene:
				key = vgene.group(1) + "_" + jgene.group(1)
				key = re.sub("[()/]","",key) #so /OR or (II) genes don't screw up the file system
				if key not in vj_partition:
					temp = "%s/%s.fa"%(prj_tree.lineage, key)
					vj_partition[key] = { 'group':key, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

				vj_partition[key]['count'] += 1
				vj_partition[key]['ids'].append(sequence.id)
				cdr3_info[sequence.id] = { 'genes' : key, 'cdr3_seq' : sequence.seq }

				#make sizes available to vsearch
				sequence.id += ";size=%d" % seqSize[sequence.id] #do this even if there's no label
										 #so I don't need to divide the cases for vsearch
				#and write
				SeqIO.write([sequence], vj_partition[key]['handle'], 'fasta')
			else:
				print("Couldn't find V and J genes for %s %s, skipping..." % (sequence.id, sequence.description))


	global natives
	natives = dict()
	natCellList = []
	if arguments['--natives'] is not None:
		natives = load_fastas(arguments['--natives'])
		for n, s in natives.items():
			if arguments['-v'] is not None:
				key = arguments['-v'] + "_" + arguments['-j']
			else:
				v_gene_pat = re.compile("(?:v_call|V_gene)=([^*\s]+)")
				j_gene_pat = re.compile("(?:j_call|J_gene)=([^*\s]+)")
				vgene = re.search(v_gene_pat, s.description)
				jgene = re.search(j_gene_pat, s.description)
				if vgene and jgene:
					key = vgene.group(1) + "_" + jgene.group(1)
				else:
					sys.exit("Can't find V and J gene annotations for native sequence %s. Please specify using the -v and -j parameters." % n)

			key = re.sub("[()/]","",key) #wouldn't expect this to be relevant for natives, but just in case...

			if key not in vj_partition:
				temp = "%s/%s.fa"%(prj_tree.lineage, key)
				vj_partition[key] = { 'group':key, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

			if arguments['--singlecell']:
				natCellID = re.search( "cell_id=(\S+)", s.description )
				if natCellID:
					cell_dict[ natCellID.group(1) ].append(n)
					natCellList.append(natCellID.group(1))
				else:
					sys.exit( f"Couldn't find `cell_id=*` annotation for native sequence {n}. Please add this information to the fasta def line." )


			seqSize[ n ] = 0
			s.id += ";size=1"
			vj_partition[key]['count'] += 1
			vj_partition[key]['ids'].append( n )
			cdr3_info[ n ] = { 'genes': key, 'cdr3_seq' : s.seq }
			SeqIO.write([ s ], vj_partition[key]['handle'], 'fasta')

	#close the file handles and delete the reference, so dict can be pickled for multithreading
	for cluster in vj_partition:
		vj_partition[cluster]['handle'].close()
		del vj_partition[cluster]['handle']

	#now go through and cluster each V/J grouping
	clusterLookup = dict()
	centroidData = dict()
	clusterSizes = Counter()
	if arguments['-t'] > 1:
		pool = Pool(arguments['-t'])
		blob = pool.map( processClusters, iterator_slice(vj_partition.values(), 25) ) #number per slice needs optimization
		pool.close()
		pool.join()

		for d in blob:
			clusterLookup.update(d['cl'])
			centroidData.update(d['cd'])
			clusterSizes.update(d['cs'])
	else:
		#don't thread
		d = processClusters( (0, vj_partition.values()) )
		clusterLookup.update(d['cl'])
		centroidData.update(d['cd'])
		clusterSizes.update(d['cs'])

	if arguments['--singlecell']:
		clusterLookup,centroidData,clusterSizes = jointClonality(clusterLookup, cell_dict, cdr3_info, natCellList)

		#now process all clusters and do tabular output
		with open( "%s/%s_lineages.txt" % (prj_tree.tables, prj_name), "w" ) as handle:
			writer = csv.writer(handle, delimiter=sep)
			writer.writerow([ "clone_id", "clone_count", "mAb_count", "junctions_nt", "junctions_aa", "included_mAbs" ])
			for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
				if size < 0:
					break
				centroidData[centroid]['rank'] = "%05d" % (rank+1)
				writer.writerow([ centroidData[centroid]['rank'], size, len(centroidData[centroid]['nats']), centroidData[centroid]['nt'],
					  				centroidData[centroid]['aa'], ",".join(sorted(centroidData[centroid]['nats'])) ])

		#update the cell_stats table
		with open("updateCellStats.tsv", 'w') as outfh:
			writer = csv.writer(outfh, delimiter="\t")
			with open( f"{prj_tree.tables}/{prj_name}_cell_stats.tsv", 'r') as infh:
				reader = csv.reader(infh, delimiter="\t")
				header = next(reader)
				hasClone = True
				if 'clone' not in header:
					hasClone = False
					header.insert(2, 'clone')
				writer.writerow(header)
				for row in reader:
					clone_id = ""
					if row[0] in clusterLookup:
						clone_id = centroidData[ clusterLookup[row[0]] ]['rank']
					if hasClone:
						row[2] = clone_id
					else:
						row.insert( 2, clone_id )
					writer.writerow( row )

		os.rename( "updateCellStats.tsv", f"{prj_tree.tables}/{prj_name}_cell_stats.tsv" )


	else:

		#regular bulk sequencing

		#now process all clusters and do tabular output
		with open( "%s/%s_lineages.txt" % (prj_tree.tables, prj_name), "w" ) as handle:
			writer = csv.writer(handle, delimiter=sep)
			writer.writerow([ "clone_id", "sequence_id", "v_call", "j_call", "junction_length_aa",
					  "junction_aa", "clone_count", "mAb_count", "included_mAbs" ])
			for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
				centroidData[centroid]['rank'] = "%05d" % (rank+1)
				writer.writerow([ "%05d"%(rank+1), centroid, centroidData[centroid]['vgene'], centroidData[centroid]['jgene'],
						  int(len(cdr3_info[centroid]['cdr3_seq'])/3), cdr3_info[centroid]['cdr3_seq'].translate(), size, len(centroidData[centroid]['nats']), ",".join(centroidData[centroid]['nats']) ])

		#do sequence output
		notationFile = re.sub( "\.f.+", "_lineageNotations.fa", arguments['--full'] )
		repFile	     = re.sub( "\.f.+", "_lineageRepresentatives.fa", arguments['--full'] )

		rep_seqs = []
		with open( notationFile, "w" ) as handle:
			for read in generate_read_fasta(arguments['--full']):
				if ";" in read.id:
					read.id = read.id[0:8] #this is for raw VSearch output with size annotations
							       #shouldn't be relevant in pipeline context
				if read.id not in clusterLookup: continue
				read.description += " clone_id=%s clone_rep=%s clone_count=%d" % ( centroidData[clusterLookup[read.id]]['rank'],
														   clusterLookup[read.id], clusterSizes[clusterLookup[read.id]] )
				SeqIO.write([read],handle,"fasta")
				if read.id in centroidData:
					rep_seqs.append(read)

		with open( repFile, "w" ) as handle:
			#use a sort to put them out in order of lineage rank (ie size)
			SeqIO.write( sorted(rep_seqs, key=lambda cent: centroidData[cent.id]['rank']), handle, "fasta" )


	#do AIRR output (both bulk and single cell)
	if arguments['--singlecell'] or "output/sequences/nucleotide" in arguments['--full']:
		if arguments['--singlecell']:
			original = "%s/%s_rearrangements_single-cell.tsv"%(prj_tree.tables, prj_name)
		elif os.path.isfile("%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name)):
			original = "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name)
		else:
			sys.exit(0) #no rearrangements file, all done

		withLin = airr.derive_rearrangement( "updateRearrangements.tsv", original, fields=["clone_id", "clone_count"])
		for r in airr.read_rearrangement( original ):
			if r['sequence_id'] in clusterLookup:
				r['clone_id']	 = centroidData[ clusterLookup[ r['sequence_id'] ] ][ 'rank' ]
				r['clone_count'] = clusterSizes[ clusterLookup[ r['sequence_id'] ] ]
			elif arguments['--singlecell'] and r['cell_id'] in clusterLookup:
				r['clone_id']	 = centroidData[ clusterLookup[ r['cell_id'] ] ][ 'rank' ]
				r['clone_count'] = clusterSizes[ clusterLookup[ r['cell_id'] ] ]
			else:
				#prevent mix-and-match data if this gets run multiple times with multiple settings
				r['clone_id']	 = ""
				r['clone_count'] = ""

			withLin.write(r)
		withLin.close()
		os.rename( "updateRearrangements.tsv", original )



if __name__ == '__main__':

	arguments = docopt(__doc__)

	arguments['--id']   = int( arguments['--id'] )
	arguments['--gaps'] = int( arguments['--gaps'] )
	arguments['-t']     = int( arguments['-t'] )

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

	if arguments['--singlecell']:

		#check for the networkx package
		try:
			from networkx import Graph, find_cliques
		except ModuleNotFoundError:
			sys.exit("The networkx package is required for single cell clonal clustering.\nPlease run `pip3 install networkx --user` and then try again.")

		#check for AIRR TSV
		if not os.path.isfile(f"{prj_tree.tables}/{prj_name}_rearrangements_single-cell.tsv"):
			sys.exit(f"{prj_tree.tables}/{prj_name}_rearrangements_single-cell.tsv does not exist.\nPlease run 1.5-single_cell_statistics.py and try again.")


	#regular (old) way
	elif arguments['--full'] is not None:
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
