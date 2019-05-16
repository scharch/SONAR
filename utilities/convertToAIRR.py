#!/usr/bin/env python3

"""
convertToAIRR.py

This script updates an _all_seq_stats.txt file from older versions of SONAR to
    an AIRR-compatible _rearrangements.tsv file used by SONAR v3.

Usage: convertToAIRR.py

Created by Chaim A. Schramm 2018-10-18.

Copyright (c) 2011-2018 Vaccine Research Center, National Institutes of Health,
                   USA. All rights reserved.

"""

import csv
from collections import defaultdict
from docopt import docopt
import airr

try:
	from SONAR import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/utilities")
	sys.path.append(find_SONAR[0])
	from SONAR import *



def main():

	airrFile = airr.create_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name), fields=['vj_in_frame','stop_codon','locus','c_call','junction_length','source_file','source_id','length_raw','length_trimmed','indels','status','blast_identity','cluster_count','v_identity'])

	
	#try to vacuum up all possible raw sequences and hope it doesn't kill memory
	raw_seqs  = defaultdict( dict )
	file_list = glob.glob("*.fa") + glob.glob("*.fas") + glob.glob("*.fst") + glob.glob("*.fasta") + glob.glob("*.fna") + glob.glob("*.fq") + glob.glob("*.fastq")
	for myseq, myqual, file_name in generate_read_fasta_folder( file_list ):
		raw_seqs[file_name][myseq.seq_id] = myseq.seq


	#get trimmed sequences
	trim_seqs = load_fastas( "%s/%s_allJ.fa"%(prj_tree.nt, prj_name) )

	#get nt junctions
	junc_seqs = load_fastas( "%s/%s_allCDR3.fa"%(prj_tree.nt, prj_name) )


	#do the conversion
	with open( "%s/%s_all_seq_stats.txt"%(prj_tree.tables, prj_name), "r" ) as handle:
		oldFile = csv.reader( handle, delimiter="\t" )
		header = next(oldFile)
		for row in oldFile:
			if row[11] == "wrong_length":
				continue

			if row[1] not in raw_seqs:
				sys.stderr.write("Couldn't find raw sequence file %s, %s will be dropped from converted file.\n"%(row[1],row[0]))
				continue
			elif row[2] not in raw_seqs[row[1]]:
				sys.stderr.write("Couldn't find raw sequence %s in file %s; %s will be dropped from converted file.\n"%(row[2],row[1],row[0]))
				continue

			r = dict()

			r['sequence']		= raw_seqs[ row[1] ][ row[2] ]
			r['sequence_alignment'] = str( trim_seqs.get( row[0], SeqRecord(seq="") ).seq )
			r['junction']		= str( junc_seqs.get( row[0], SeqRecord(seq="") ).seq )
					
			r['sequence_id']	= row[0]
			r['source_file']	= row[1]
			r['source_id']		= row[2]
			r['length_raw']		= row[3]
			if not row[4] == "NA":
				r['length_trimmed']  = row[4]
			if not row[5] == "NA":
				r['v_call']	     = row[5]
			if row[6] not in ["NA", "not_found"]:
				r['d_call']	     = row[6]
			if not row[7] == "NA":
				r['j_call']	     = row[7]
			if not row[9] == "NA":
				r['indels']	     = row[9]
			if not row[10] == "NA":
				r['stop_codon']	     = row[10]
			r['status']		     = row[11]
			if not row[12] == "NA":
				r['blast_identity']  = "%.3f" % (  1 - float(re.sub("%","",row[12]))/100  )
			if not row[13] == "NA":
				r['junction_length'] = int(row[13])+6
			if not row[15] == "NA":
				r['junction_aa']     = row[15]

			if len(row)>15:
				if header[16]=="Unique":
					if row[16] == "T":
						r['status']	   = "unique"
						r['cluster_count'] = row[17]
					if len(row)>17 and not row[18]=="NA":
						r['v_identity']	   = "%.3f" % (	 1 - float(re.sub("%","",row[18]))/100	)
				elif header[16] == "V_div" and not row[16]=="NA":
					r['v_identity']		   = "%.3f" % (	 1 - float(re.sub("%","",row[16]))/100	)


			#figure out in-frame/productive
			if row[10] == "good":
				r['vj_in_frame'] = "T"
				r['productive']	 = "T"
			elif row[10] == "stop":
				r['vj_in_frame'] = "T"
				r['productive']	 = "F"
			elif row[10] == "nonproductive":
				r['vj_in_frame'] = "F"
				r['productive']	 = "F"
			elif row[10] == "indel":
				r['productive']	 = "F"

				
			#figure out locus
			if any( x in row[5] for x in ["HV", "VH", "Vh", "vh", "heavy", "Heavy", "HEAVY"] ):
				r['locus'] = "IGH"
			elif any( x in row[5] for x in ["LV", "VL", "Vl", "vl", "lambda", "Lambda", "LAMBDA"] ):
				r['locus'] = "IGL"
			elif any( x in row[5] for x in ["KV", "VK", "Vk", "vk", "kappa", "Kappa", "KAPPA"] ):
				r['locus'] = "IGK"

				
			airrFile.write(r)

	airrFile.close()
	valid = airr.validate_rearrangement( "%s/%s_rearrangements.tsv"%(prj_tree.tables, prj_name) )
	if not valid:
		sys.exit( "ERROR: something went wrong, %s/%s_rearrangements.tsv failed validation!"%(prj_tree.tables, prj_name) )

		
			
if __name__ == "__main__":

	#parse arguments	
	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	if not os.path.isfile( "%s/%s_all_seq_stats.txt"%(prj_tree.tables, prj_name) ):
		sys.exit( "No old file to convert!" )
		   
	#log command line
	logCmdLine(sys.argv)

	main()
