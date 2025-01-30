#!/usr/bin/env python3

"""
2.4-cluster_into_groups.py

This script uses CDR3 identity to group unique sequences from a given data set
      into pseudo-lineages that can help define groups of related B cells.
      Uses the AIRR-formatted rearrangements.tsv as default input. FASTA input
      is no longer accepted, but multiple input rearrangements TSVs can be used
      to add sequences from other sources, like experimentally isolated
      monoclonals.

      Sequences are first grouped by unique V and J gene assignments and then
      VSearch is used to cluster the CDR3 sequences. Alternatively, closely 
      related V genes may be treated as indistinguishable when partitioning 
      reads for clustering, using the `--geneClusters` flag. This is based
      on the approach of Luo, Yu, and Song, PLoS Comp Biol 2016. Predetermined
      clusters are provided for default databases using heirarchical
      clustering with a 4% threshold, based on an analysis of interderminate
      gene assignments in a variety of historical datasets. Custom clusters
      can be supplied using `--customClusters`.

      The default threshold of 90% identity (bulk) or 80% identity (single cell)
      with no in-dels is probably useful for most cases, but more stringent or 
      lenient criteria may sometimes be more appropriate. 

Usage: 2.4-cluster_into_groups.py [ --rearrangements TSV... --names SAMPLE... --filter all --id <90> --gaps <0> --output TSV --geneClusters --customClusters <clusters.txt> --species <human> --singlecell --preserve --master <db.xlsx> --subject A123 -t 1 ]

Options:
    --rearrangements TSV               One or more AIRR-formatted rearrangements files with the 
                                          sequences to be clustered into lineages.
                                          [default: output/tables/<project>_rearrangements.tsv]
    --names SAMPLE                     Optional short names to keep track of which of multiple 
                                          input files output rearrangements are derived from. If 
                                          specified, the number of names provided *must* match the 
                                          number of input `rearrangements` files. If no `names` are
                                          given, the full paths specified to `rearrangements` will 
                                          be used. As a special case, the usage `--names preserve` 
                                          will extract the short names from an existing 
                                          `source_repertoire` column.
    --filter all                       Filter sequences by status before calculating lineages. 
                                          Allowed values are "all" (ie all CDR3), "good", "unique", 
                                          (determined by having `centroid`==`sequence_id` --does 
                                          NOT remove singletons!) and "paired" (with `--singlecell`
                                          only).  [default: all]
    --id <90>                          Clustering threshold (%) for CDR3 sequence identity 
                                          (nucleotide). (default: 90, or 80 with `--singlecell`)
    --gaps <0>                         Maximum number of in-dels to allow between CDR3 sequences in
                                          the same group. NOTE: This is implemented in nucleotide 
                                          space using vsearch's --maxgaps parameter, which counts 
                                          gap openings, rather than gap columns, so probably 
                                          shouldn't be set to more than 1. There is no guarantee an
                                          even codon in-del in the alignment used for clustering! 
                                          Also, vsearch will still count the gap columns as 
                                          mismatches, so a CDR3 of 20AA will be counted as 95% id 
                                          to an identical-other-than-deletion 19AA CDR3. Set your
                                          threshold for --id accordingly. [default: 0]
    --output TSV                       File where the output should be saved. If not specified, 
                                          output will overwrite the first input file.
    --geneClusters                     Flag to indicate that reads should be partitioned based on 
                                          closely related V genes that may be prone to mutual 
                                          misassignment, instead of using exact matches of assigned
                                          V and J genes. Using this option will turn off matching 
                                          on J genes entirely. Predetermined clusters for supported
                                          species are in SONAR/sample_data/functionalClusters. 
                                          Otherwise, use `--customClusters` to specify.
    --customClusters <clusters.txt>    Option specifiying a file with gene clusters for non-default
                                          gene databases. Format is a tab-delimited file with one
                                          line for each allele in the V gene database (first
                                          column), with the cluster name/id indicated in the second
                                          column. Only relevant if `--geneClusters` is used.
                                          Takes priority over default clusters (and `--species`).
    --species <human>                  Option to indicate the (supported) species to use for 
                                          geneClusters if the clonality analysis is being done in a
                                          different directory than the original annotation (ie. 
                                          work/internal/gene_locus.txt is missing).
    --singlecell                       A flag to indicate single cell data - heavy and light chain 
                                          data will be used jointly to define clones. `cell_id` 
                                          column must be present in all input rearrangements files.
                                          Note that if the `cell_status` column is present, 
                                          suspected multiplets will be filtered out.
    --preserve                         A flag to preserve the original clone IDs when adding new 
                                          sequences to previously processed data. Use with caution,
                                          as it will silently join and/or split clones even if the 
                                          actual underlying clustering came out differently. Only 
                                          operates on the first TSV if multiple are provided.
    --master <db.xlsx>                 An Excel file with a database of previously seen single cell 
                                          sequences. Extracts those corresponding to the current
                                          subject (`--subject` is required) and adds them to the
                                          lineage clustering to maintain clone IDs. Only data from
                                          the `--rearrangements` file will be included in the output.
    --subject A123                     Subject ID to extract from master database> Required if
                                          `--master` is used, ignored otherwise.
    -t 1                               Number of threads used [default: 1]


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
Switched over to AIRR TSV input only by CAS 2020-05-19.
Added --preserve option by CAS 2020-07-02.
Updated filters to new syntax by CAS 2020-07-02.
Added geneClusters option by CAS on 2020-07-16.
Added short names option by CAS on 2020-07-16.
Added numRepertoires column to lineage output by CAS, 2020-07-16.
Fixed short names by CAS in 2020-07-28.
Added ability to get short names from the input file to allow adding new data
                         multiple times by CAS 2020-10-22.
Fixed new clone numbers when using --preserve by CAS 2020-10-22.
Fixed representative CDR3s for single cell clustering by CAS 2020-10-22. 
Added support for --species and --customClusters to make --geneClusters
                         more flexible/useful by CA Schramm 2022-05-16.
Changed default id threshold for single cells to 80% by CAS 2022-07-14.
Added clean up to end of run by CA Schramm 2022-07-14.
Lineage table will print actual V genes instead of OTU centroid when in
                         geneClusters mode by CA Schramm 2022-07-14.
Added option to pull old sequences directly from master db to minimize
                         errors by CA Schramm 2025-01-09.
Switched to handling of cell_stats files via pandas to keep columns from
                         getting messed up by CA Schramm 2025-01-10.

Copyright (c) 2011-2025 Columbia University and Vaccine Research Center, National
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
import pandas as pd

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

	clusterLookup = dict()
	centroidData = dict()
	clusterSizes = Counter()
	countsByInput = defaultdict( Counter )

	for cluster in chunk:

		#save a bit of time for obvious singletons
		if cluster['count'] == 1:
			single = cluster['ids'][0]
			myGenes = cluster['group'].split("_")
			clusterLookup[ single ] = single
			if arguments['--geneClusters']:
				centroidData[ single ] = dict( vgene = myGenes[0], jgene = "" )
			else:
				centroidData[ single ] = dict( vgene = myGenes[0], jgene = myGenes[1] )
			clusterSizes[ single ] = 1#seqSize[ single ]
			origin = re.search("===(.+)$", single).groups()[0]
			countsByInput[ single ][ origin ] += 1
			continue

		#cluster with vsearch
		subprocess.call([vsearch, "-cluster_size", cluster['file'],
				 "-id", str(arguments['--id']/100.0),
				 "-maxgaps", str(arguments['--gaps']),
				 "-sizein", "-uc", "%s/%s.uc"%(prj_tree.lineage, cluster['group']),
				 "-minseqlength", "15", #lets us capture CDR3s down to 3 aa
				 "-wordlength", "3", #don't require long homologous blocks
				 "-minwordmatches", "1", #turn sensitivity all the way up, rely on percent id for specificity
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

				origin = re.search("===(.+)$", hit).groups()[0]

				if row[0] == "S":
					if arguments['--geneClusters']:
						centroidData[ hit ] = dict( vgene = myGenes[0], jgene = "" )
					else:
						centroidData[ hit ] = dict( vgene = myGenes[0], jgene = myGenes[1] )
					clusterSizes[ hit ] = 1#seqSize[ hit ]
					clusterLookup[ hit ] = hit
					countsByInput[ hit ][ origin ] += 1
				elif row[0] == "H":
					clusterLookup[ hit ] = cent
					clusterSizes[ cent ] += 1#seqSize[ hit ]
					countsByInput[ cent ][ origin ] += 1
				else:
					break #skip "C" lines

	return { 'cl':clusterLookup, 'cd':centroidData, 'cs':clusterSizes, 'ci': countsByInput }



def jointClonality(clusters, cells, cdr3Info):

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
	cellClones = dict()
	for cliq in find_cliques(cloneGraph):
		cellClones[ frozenset(cliq) ] = []

	#Go back through the list of cells and assign them to clones
	assignments = dict()
	countsByInput = defaultdict( Counter )
	ambiguous = 0
	for c in cells:
		origin = re.search("===(.+)$", c).groups()[0]
		possibleClones = []
		for cc in cellClones:
			if frozenset([clusters[s] for s in cells[c]]) <= cc:
				possibleClones.append(cc)
		if len(possibleClones) == 1:
				assignments[c] = possibleClones[0]
				cellClones[ possibleClones[0] ].append(c)
				countsByInput[ possibleClones[0] ][ origin ] += 1
		else:
			#no matches found, probably ambiguous
			# but check if its a singleton that never appeared with any other chain
			if len(possibleClones) == 0 and len(cells[c]) == 1:
				cloneName = frozenset([clusters[s] for s in cells[c]])
				assignments[c] = cloneName
				cellClones[ cloneName ] = [c]
				countsByInput[ cloneName ][ origin ] += 1
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
			#The centroid of the original clustering (ie `chain`) is not guaranteed to be
			#    in the final `cellClone`, especially at low identity thresholds. Rather
			#    than trying to recluster to get a new centroid, just grab the first one:
			found = False
			for myCell in cellClones[cc]:
				for myChain in cells[myCell]:
					if cdr3Info[myChain]['genes'] == cdr3Info[chain]['genes']:
						nt.append( f"{cdr3Info[myChain]['genes']}:{cdr3Info[myChain]['cdr3_seq']}" )
						aa.append( f"{cdr3Info[myChain]['genes']}:{cdr3Info[myChain]['cdr3_seq'].translate()}" )
						found = True
						break
				if found:
					break
		cloneInfo[cc] = { 'aa':",".join(sorted([x.upper() for x in aa])), 'nt':",".join(sorted([x.upper() for x in nt])) }

	#turn cellClones into a Counter for compatibility with code for bulk sequencing
	cellCloneCounter = Counter()
	for cc in cellClones:
		cellCloneCounter[cc] = len(cellClones[cc])
	return assignments, cloneInfo, cellCloneCounter, countsByInput



def main():

	#decide which categories of reads we are going to include
	#first for cells, then for individual rearrangements
	filter_rules = []
	if arguments['--singlecell']:
		filter_rules.append( 're.search("^(?:(?!multi|none).)*$", r["cell_status"])' )
		if arguments['--filter'] == "paired":
			filter_rules = [ "r['cell_status'] == 'canonical_pair'" ]
	if arguments['--filter'] == "all":
		filter_rules.append( "r['junction'] != ''" )
	elif arguments['--filter'] == "good":
		filter_rules.append( "r['status'] == 'good'")
	elif arguments['--filter'] == "unique":
		filter_rules.append( "r['centroid'] == r['sequence_id']" )

	vj_partition = dict()
	cdr3_info = dict()
	seqSize = Counter()
	oldClones = dict()
	sourceList = list()
	cell_dict = defaultdict(list)

	#extract sequences from master database
	if arguments['--master'] is not None:
		wb = openpyxl.load_workbook( arguments['--master'], read_only=True )
		ws = wb.active
		ws.reset_dimensions()
		data = ws.values
		header = next(data)
		suCol = header.index("SubjectID")
		lnCol = header.index("LineageNum")
		abCol = header.index("AntibodyNum")
		vhCol = header.index("VH_call")
		jhCol = header.index("JH_call")
		vlCol = header.index("VKL_call")
		jlCol = header.index("JKL_call")
		h3Col = header.index("CDRH3_nt")
		aaCol = header.index("VJ_aa")
		ntCol = header.index("VJ_nt")
		l3Col = header.index("CDRL3_AA")
		dbSeqs = 0
		for row in ws.values:
			if row[suCol] == arguments['--subject']:
				dbSeqs += 1
				cell = f"{row[suCol]}-{row[lnCol]}.{row[abCol]}===masterDB"
				hSeq = f"H-{cell}"
				lSeq = f"L-{cell}"
				cdrl3_seq = row[ntCol][ row[aaCol].index(row[l3Col])*3 : (row[aaCol].index(row[l3Col])+len(row[l3Col]))*3 ]

				cell_dict[ cell ] += [hSeq, lSeq]
				oldClones[ cell ] = row[ lnCol ]

				keyH = row[vhCol].split("*")[0] + "_" + row[jhCol].split("*")[0]
				cdr3_info[ hSeq ] = { 'genes' : keyH, 'cdr3_seq' : Seq.Seq(row[h3Col]) }
				keyL = row[vlCol].split("*")[0] + "_" + row[jlCol].split("*")[0]
				cdr3_info[ lSeq ] = { 'genes' : keyL, 'cdr3_seq' : Seq.Seq(cdrl3_seq) }

				if arguments['--geneClusters']:
					keyH = geneClusters.get( row[vhCol].split(",")[0], row[vhCol].split("*")[0] )
					cdr3_info[ hSeq ] = { 'genes' : row[vhCol].split("*")[0], 'cdr3_seq' : Seq.Seq(row[h3Col]) }
					keyL = geneClusters.get( row[vlCol].split(",")[0], row[vlCol].split("*")[0] )
					cdr3_info[ lSeq ] = { 'genes' : row[vlCol].split("*")[0], 'cdr3_seq' : Seq.Seq(cdrl3_seq) }

				if keyH not in vj_partition:
					temp = "%s/%s.fa"%(prj_tree.lineage, keyH)
					vj_partition[keyH] = { 'group':keyH, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }
				if keyL not in vj_partition:
					temp = "%s/%s.fa"%(prj_tree.lineage, keyL)
					vj_partition[keyL] = { 'group':keyL, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

				vj_partition[keyH]['count'] += 1
				vj_partition[keyH]['ids'].append(hSeq)
				vj_partition[keyL]['count'] += 1
				vj_partition[keyL]['ids'].append(lSeq)

				#create a sequence object and write
				tempSeqH = SeqRecord( id=hSeq, seq=Seq.Seq(re.sub("[-.+]","",row[h3Col])) )
				tempSeqH.id += ";size=1"
				SeqIO.write([tempSeqH], vj_partition[keyH]['handle'], 'fasta')

				tempSeqL = SeqRecord( id=lSeq, seq=Seq.Seq(re.sub("[-.+]","",cdrl3_seq)) )
				tempSeqL.id += ";size=1"
				SeqIO.write([tempSeqL], vj_partition[keyL]['handle'], 'fasta')

		wb.close()
		if dbSeqs == 0:
			print( f"Warning: No rows for subject {arguments['--subject']} found in master database {arguments['--master']}" )

	#open the rearrangement file(s) and parse into groups with same V/J
	for index, inFile in enumerate(arguments['--rearrangements']):

		#open the file and check that we have all the required fields
		reader = airr.read_rearrangement(inFile)
		if arguments['--singlecell'] and not "cell_id" in reader.fields:
			sys.exit( f"`cell_id` column not found in {inFile}, cannot do single-cell lineage analysis" )
		if arguments['--filter'] == "unique" and not "centroid" in reader.external_fields:
			sys.exit( f"`centroid` column not found in {inFile}, cannot filter for unique sequences" )
		if arguments['--preserve'] and arguments['--master'] is None and index==0 and not "clone_id" in reader.fields:
			print("Can't find existing clone_ids, `--preserve` will be ignored...", file=sys.stderr)
		if len(arguments['--names']) > 0 and arguments['--names'][index] == "preserve" and not "source_repertoire" in reader.fields:
			sys.exit("Can't find existing source_repertoire, please use a `--name` other than 'preserve'." )

		#now iterate through the rearrangements
		for r in filterAirrTsv(reader, filter_rules):

			#short CDR3s can't be clustered and probably indicate a bad sequence anyway
			if r['junction_length'] < 15:
				continue

			#skip sequences with missing gene assignments
			if r['v_call'] == "" or r['j_call'] == "":
				continue

			#uniquify the sequence id in case of multiple file inputs
			suffix = index
			if len(arguments['--names']) > 0:
				if arguments['--names'][index] == "preserve":
					suffix = r[ 'source_repertoire' ]
				else:
					suffix = arguments['--names'][index]
			if suffix not in sourceList:
				sourceList.append(suffix)
			r[ 'sequence_id' ] += f"==={suffix}"

			if arguments['--singlecell']:
				#for network analysis
				r[ 'cell_id' ] += f"==={suffix}"
				cell_dict[ r['cell_id'] ].append(r['sequence_id'])

			if arguments['--preserve'] and arguments['--master'] is None and index==0 and 'clone_id' in r and r['clone_id']!="":
				if arguments['--singlecell']:
					oldClones[ r['cell_id'] ] = r['clone_id']
				else:
					oldClones[ r['sequence_id'] ] = r['clone_id']

			#get size
			seqSize[ r['sequence_id'] ] = 1
			if 'consensus_count' in r and r['consensus_count'] is not None and not r['consensus_count'] == "":
				seqSize[ r['sequence_id'] ] = r['consensus_count']
			elif 'duplicate_count' in r and r['duplicate_count'] is not None and not r['duplicate_count'] == "":
				seqSize[ r['sequence_id'] ] = r['duplicate_count']

			#get gene assignments
			key = r['v_call'].split("*")[0] + "_" + r['j_call'].split("*")[0]
			cdr3_info[ r['sequence_id'] ] = { 'genes' : key, 'cdr3_seq' : Seq.Seq(r['junction']) }

			if arguments['--geneClusters']:
				key = geneClusters.get( r['v_call'].split(",")[0], r['v_call'].split("*")[0] )
				cdr3_info[ r['sequence_id'] ] = { 'genes' : r['v_call'].split("*")[0], 'cdr3_seq' : Seq.Seq(r['junction']) }

			if key not in vj_partition:
				temp = "%s/%s.fa"%(prj_tree.lineage, key)
				vj_partition[key] = { 'group':key, 'handle':open(temp, "w"), 'file':temp, 'count':0, 'ids':[] }

			vj_partition[key]['count'] += 1
			vj_partition[key]['ids'].append(r['sequence_id'])

			#create a sequence object
			tempSeq = SeqRecord( id=r['sequence_id'], seq=Seq.Seq(re.sub("[-.+]","",r['junction'])) )
			tempSeq.id += ";size=%d" % seqSize[ r['sequence_id'] ] #do this even if there's no label
									                               #so I don't need to divide the cases for vsearch
			#and write
			SeqIO.write([tempSeq], vj_partition[key]['handle'], 'fasta')


	#close the file handles and delete the reference, so dict can be pickled for multithreading
	for cluster in vj_partition:
		vj_partition[cluster]['handle'].close()
		del vj_partition[cluster]['handle']

	#now go through and cluster each V/J grouping
	clusterLookup = dict()
	centroidData = dict()
	clusterSizes = Counter()
	countsByInput = defaultdict( Counter )
	if arguments['-t'] > 1:
		pool = Pool(arguments['-t'])
		blob = pool.map( processClusters, iterator_slice(vj_partition.values(), 25) ) #number per slice needs optimization
		pool.close()
		pool.join()

		for d in blob:
			clusterLookup.update(d['cl'])
			centroidData.update(d['cd'])
			clusterSizes.update(d['cs'])
			countsByInput.update(d['ci'])
	else:
		#don't thread
		d = processClusters( (0, vj_partition.values()) )
		clusterLookup.update(d['cl'])
		centroidData.update(d['cd'])
		clusterSizes.update(d['cs'])
		countsByInput.update(d['ci'])

	#make some output file names
	lineageFile  = re.sub("(_rearrangements.*)?\.tsv", "_lineages.txt", arguments['--output'])
	cellStatFile = re.sub("(_rearrangements.*)?\.tsv", "_cell_stats.tsv", arguments['--output'])
	#make sure we don't accidentally overwrite anything
	if lineageFile == arguments['--output']:
		lineageFile  = arguments['--output'] + "_lineages.txt"
		cellStatFile = arguments['--output'] + "_cell_stats.tsv"


	#do joint clonality for single cells
	if arguments['--singlecell']:
		clusterLookup,centroidData,clusterSizes,countsByInput = jointClonality(clusterLookup, cell_dict, cdr3_info)

		#now process all clusters and do tabular output
		with open( lineageFile, "w" ) as handle:
			writer = csv.writer(handle, delimiter=sep, dialect='unix', quoting=csv.QUOTE_NONE)

			header = [ "clone_id", "clone_count", "junctions_nt", "junctions_aa" ]
			if len(arguments['--rearrangements']) > 1:
				header += ["source_count", "num_sources"]
			writer.writerow(header)

			#if we are reclustering with the --preserve option, we need to figure out
			#    the numbering
			currentMaxCloneNum = len(set(oldClones.values()))
			if len(oldClones.keys()) > 0:
				currentMaxCloneNum = max( map(int, oldClones.values()) )

			for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
				if size == 0:
					break

				if arguments['--preserve']:
					#get cells with this centroid to look up old clone_id
					oldCells = [k for k,v in clusterLookup.items() if v == centroid and k in oldClones]

					if len(oldCells)>0:
						centroidData[centroid]['rank'] = oldClones[ oldCells[0] ]
					else:
						currentMaxCloneNum += 1
						centroidData[centroid]['rank'] = "%05d" % currentMaxCloneNum
				else:
					centroidData[centroid]['rank'] = "%05d" % (rank+1)

				dataToWrite = [ centroidData[centroid]['rank'], size, centroidData[centroid]['nt'],
					  				centroidData[centroid]['aa'] ]
		
				#find how many members are from each source file
				if len(sourceList) > 1:
					breakdown = [ countsByInput[centroid][suffix] for suffix in sourceList ]
					dataToWrite += [ ":".join([ str(b) for b in breakdown]), sum([1 if b>0 else 0 for b in breakdown]) ]

				writer.writerow( dataToWrite )

		#update the cell_stats table
		#for each input try to guess the matching cell_stats file
		allDFs = []
		for ind in range(len(arguments['--rearrangements'])):
			cell_stats = re.sub("_rearrangements.*\.tsv", "_cell_stats.tsv", arguments['--rearrangements'][ind])

			airrFile = arguments['--rearrangements'][ind]
			if len(arguments['--names']) > 0:
				airrFile = arguments['--names'][ind]

			#check if it exists
			if not os.path.isfile(cell_stats):
				print( f"Warning: Cannot find cell stats file for {arguments['--rearrangements'][ind]} (tried {cell_stats}).\nCells from this repertoire will not be included in output to {cellStatFile}\n", file=sys.stderr)
				continue

			df = pd.read_csv(cell_stats, sep="\t", dtype={'clone':'str'})
			if 'clone' not in df.columns:
				df.insert(2,"clone","")
			if 'source' not in df.columns:
				df.insert(3,"source","")

			for index,row in df.iterrows():
				clone_id = ""
				suffix = ind
				if len(arguments['--names']) > 0:
					if arguments['--names'][ind] == "preserve":
						suffix = row['clone']
					else:
						suffix = arguments['--names'][ind]
				try:
					unique_cell = row['cell']+f"==={suffix}"
				except TypeError:
					sys.exit( f"Please check alignment of columns and make sure there are sufficient column titles in the header for {cell_stats}")
				if unique_cell in clusterLookup:
					clone_id = centroidData[ clusterLookup[ unique_cell ] ]['rank']
				df.at[index, 'clone'] = clone_id

				if len(arguments['--rearrangements']) > 1 or len(arguments['--names']) > 0:
					if airrFile != "preserve":
							df.at[index, 'source'] = airrFile

			allDFs.append(df)

		bigdf = pd.concat( allDFs, ignore_index=True )
		bigdf.to_csv(cellStatFile, index=False, sep="\t")

	else:

		#regular bulk sequencing

		#now process all clusters and do tabular output
		with open( lineageFile, "w" ) as handle:
			writer = csv.writer(handle, delimiter=sep, dialect='unix', quoting=csv.QUOTE_NONE)
			header = [ "clone_id", "centroid", "v_call", "j_call", "junction_length_aa",
					  "junction_aa", "clone_count" ]
			if len(arguments['--rearrangements']) > 1:
				header += ['source_count', 'num_sources']
			writer.writerow(header)

			for rank, (centroid, size) in enumerate(clusterSizes.most_common()):
				if arguments['--preserve']:
					if centroid in oldClones:
						centroidData[centroid]['rank'] = oldClones[centroid]
					else:
						centroidData[centroid]['rank'] = "%05d" % (len(oldClones.keys()) + rank)
				else:
					centroidData[centroid]['rank'] = "%05d" % (rank+1)

				dataToWrite = [ "%05d"%(rank+1), re.search("(.+)===.+$", centroid).groups()[0], centroidData[centroid]['vgene'], centroidData[centroid]['jgene'],
						  int(len(cdr3_info[centroid]['cdr3_seq'])/3), cdr3_info[centroid]['cdr3_seq'].translate(), size ]
				#find how many members are from each source file
				if len(arguments['--rearrangements']) > 0:
					breakdown = [ countsByInput[centroid][suffix] for suffix in sourceList ]
					dataToWrite += [ ":".join([ str(b) for b in breakdown]), sum([1 if b>0 else 0 for b in breakdown]) ]
				writer.writerow(dataToWrite)

	#do AIRR output (both bulk and single cell)
	#use a temp file to avoid problems trying to overwrite an input file
	extra_fields = ["clone_id", "clone_count"]
	if len(arguments['--rearrangements']) > 1 or len(arguments['--names']) > 0:
		extra_fields += ['source_repertoire']

	withLin = airr.derive_rearrangement( "temp.tsv", arguments['--rearrangements'][0], fields=extra_fields)
	for index, inFile in enumerate(arguments['--rearrangements']):
		for r in airr.read_rearrangement( inFile ):

			#get the uniquified ids
			suffix = index
			if len(arguments['--names']) > 0:
				if arguments['--names'][index] == "preserve":
					suffix = r[ 'source_repertoire' ]
				else:
					suffix = arguments['--names'][index]

			unique_seq  = r['sequence_id']     + f"==={suffix}"
			unique_cell = r.get('cell_id', '') + f"==={suffix}"

			if unique_seq in clusterLookup:
				r['clone_id']	 = centroidData[ clusterLookup[ unique_seq ] ][ 'rank' ]
				r['clone_count'] = clusterSizes[ clusterLookup[ unique_seq ] ]
			elif arguments['--singlecell'] and unique_cell in clusterLookup:
				r['clone_id']	 = centroidData[ clusterLookup[ unique_cell ] ][ 'rank' ]
				r['clone_count'] = clusterSizes[ clusterLookup[ unique_cell ] ]
			else:
				#prevent mix-and-match data if this gets run multiple times with multiple settings
				r['clone_id']	 = ""
				r['clone_count'] = ""

			#add source repertoire if relevant
			if len(arguments['--rearrangements']) > 1 or len(arguments['--names']) > 0:
				if len(arguments['--names']) > 0:
					if arguments['--names'][ index ] != "preserve":
						r['source_repertoire'] = arguments['--names'][ index ]
				else:
					r['source_repertoire'] = inFile

			withLin.write(r)
	withLin.close()

	#Put the output TSV in the desired destination
	os.rename( "temp.tsv", arguments['--output'] )

	to_clean = glob.glob("%s/*fa"%prj_tree.lineage) + glob.glob("%s/*fa"%prj_tree.lineage)
	if len(to_clean) > 0:
		print("Cleaning up...",file=sys.stderr)
		for f in to_clean:
			try:
				os.remove(f)
			except:
				pass


if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	arguments['--rearrangements'][0] = re.sub("<project>", prj_name, arguments['--rearrangements'][0])
	for airrTsv in arguments['--rearrangements']:
		if not os.path.isfile(airrTsv):
			sys.exit(f"Cannot find rearrangements file {airrTsv}")
		elif not airr.validate_rearrangement(airrTsv):
			sys.exit(f"File {airrTsv} is not in valid AIRR format.")

	if len(arguments['--names']) > 0 and len(arguments['--names']) != len(arguments['--rearrangements']):
		sys.exit("Error: number of `names` must match number of `rearrangements`.")

	if arguments['--output'] is None:
		arguments['--output'] = arguments['--rearrangements'][0]

	if arguments['--id'] is None:
		if arguments['--singlecell']:
			arguments['--id'] = 80
		else:
			arguments['--id'] = 90
	else:
		arguments['--id']   = int( arguments['--id'] )
	arguments['--gaps'] = int( arguments['--gaps'] )
	arguments['-t']     = int( arguments['-t'] )

	if not arguments['--filter'] in ["all", "good", "unique", "paired"]:
		sys.exit("Allowed values for `--filter` are 'all', 'good', 'unique', and 'paired' only.")

	geneClusters = dict()
	if arguments['--geneClusters']:

		if arguments['--customClusters'] is not None:
			if os.path.exists(arguments['--customClusters']):
				with open(arguments['--customClusters'], 'r') as database:
					reader = csv.reader(database, delimiter="\t")
					for row in reader:
						geneClusters[ row[0] ] = row[1]
			else:
				sys.exit(f"Can't find cluster file {arguments['--customClusters']}")

		elif arguments['--species'] is not None:
			if arguments['--species'] not in SUPPORTED_SPECIES:
				sys.exit( "Error: `--species` must be one of: " + ",".join(SUPPORTED_SPECIES.keys()) )
			else:
				with open(f"{SCRIPT_FOLDER}/sample_data/functionalClusters/{arguments['--species']}_clusterLookup.txt", 'r') as database:
					reader = csv.reader(database, delimiter="\t")
					for row in reader:
						geneClusters[ row[0] ] = row[1]

		else:
			if os.path.exists(f"{prj_tree.internal}/gene_locus.txt"):
				with open(f"{prj_tree.internal}/gene_locus.txt", 'r') as check:
					species = next(check).strip()
					if species in SUPPORTED_SPECIES:
						with open(f"{SCRIPT_FOLDER}/sample_data/functionalClusters/{species}_clusterLookup.txt", 'r') as database:
							reader = csv.reader(database, delimiter="\t")
							for row in reader:
								geneClusters[ row[0] ] = row[1]
					else:
						sys.exit("Please use --customClusters for external gene databases.")
			else:
				sys.exit(f"{prj_tree.internal}/gene_locus.txt not found, please use --species or --customClusters")


	if arguments['--singlecell']:
		#check for the networkx package
		try:
			from networkx import Graph, find_cliques
		except ModuleNotFoundError:
			sys.exit("The networkx package is required for single cell clonal clustering.\nPlease run `pip3 install networkx --user` and then try again.")

	if arguments['--master'] is not None:
		if not os.path.isfile(arguments['--master']):
			sys.exit(f"Cannot find database file {arguments['--master']}")
		if not arguments['--singlecell']:
			sys.exit( "Master database input is only compatible with single cell mode" )
		if arguments['--subject'] is None:
			sys.exit( "Error: must specify subject ID to be extracted from master database")
		if arguments['--preserve']:
			print( "Warning: only cloneIDs in the master database will be preserved", file=sys.stderr)
		else:
			arguments['--preserve'] = True

		#check for the openpyxl package
		try:
			import openpyxl
		except ModuleNotFoundError:
			sys.exit("The openpyxl package is required to read a Master Database.\nPlease run `pip3 install openpyxl --user` and then try again.")

	#log command line
	logCmdLine(sys.argv)


	#create necessary directories if they don't exist
	os.makedirs(prj_tree.logs, exist_ok=True)
	os.makedirs(prj_tree.lineage, exist_ok=True)
	os.makedirs(prj_tree.tables, exist_ok=True)


	main()
