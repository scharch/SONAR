#!/usr/bin/env python3

#clean up function
def cleanup():
	shutil.rmtree( "output", ignore_errors=True )
	shutil.rmtree( "work", ignore_errors=True)
	try:
		os.remove("f0_merged.fq")
		os.remove("derepAllRawSeqs.uc")
		os.remove("lineage.fa")
	except:
		pass


#test imports (it's a bit silly to include builtins in this list...)
import importlib, sys
for checkModule in ['airr', 'atexit', 'Bio', 'collections', 'colorsys', 'csv', 'datetime', 'docopt', 'ete3', 'fileinput', 'functools', 'fuzzywuzzy', 'glob', 'gzip', 'io', 'itertools', 'Levenshtein', 'math', 'multiprocessing', 'numpy', 'os', 'pandas', 'pickle', 'PyQt4.QtGui', 'random', 're', 'shutil', 'statistics', 'string', 'subprocess', 'sys', 'time', 'traceback']:

	try:
		my_module = importlib.import_module(checkModule)
	except ImportError:
		sys.exit(f"module {checkModule} not found")

#check Biopython submodules (not sure this is necessary)
for submodule in ['Align', 'AlignIO', 'Alphabet', 'Data.CodonTable', 'Phylo', 'Seq', 'SeqIO', 'SeqRecord']:
	try:
		my_module = importlib.import_module(f"Bio.{submodule}", package="Bio")
	except ImportError:
		sys.exit(f"module Bio.{submodule} not found")

#check for master script
import shutil
if not shutil.which('sonar'):
	print("Master script not found in PATH, programs will be invoked directly", file=sys.stderr)


#test SONAR scripts
import os, subprocess
SONARDIR = os.path.abspath(sys.argv[0]).split("SONAR/tests")[0]
os.chdir(f"{SONARDIR}/SONAR/tests")
for command in [ [f"{SONARDIR}/SONAR/annotate/1.0-preprocess.py", "--input", "subsample_r1.fq.gz", "--reverse", "subsample_r2.fq.gz"],
                 [f"{SONARDIR}/SONAR/annotate/1.1-blast_V.py", "--fasta", "f0_merged.fq", "--derep", "--npf", "2000", "--threads", "2"],
				 [f"{SONARDIR}/SONAR/annotate/1.2-blast_J.py"],
				 [f"{SONARDIR}/SONAR/annotate/1.3-finalize_assignments.py"],
				 [f"{SONARDIR}/SONAR/annotate/1.4-cluster_sequences.py"],
				 [f"{SONARDIR}/SONAR/lineage/2.1-calculate_id-div.py"],
				 [f"{SONARDIR}/SONAR/plotting/4.1-setup_plots.pl", "--statistic", "div"],
				 [f"{SONARDIR}/SONAR/lineage/2.4-cluster_into_groups.py"],
				 [f"{SONARDIR}/SONAR/utilities/getReadsByAnnotation.py", "-f", "output/sequences/nucleotide/tests_goodVJ_unique_lineageNotations.fa", "-a", "clone_id=000(01|07|08)", "-o", "lineage.fa"],
				 [f"{SONARDIR}/SONAR/phylogeny/3.2-run_IgPhyML.py", "-v", "IGHV4-39*01", "--seqs", "lineage.fa", "--quick", "--seed", "321325749"],
				 [f"{SONARDIR}/SONAR/utilities/flipTree.pl", "output/tests_igphyml.tree", "output/tests_igphyml.flipped.tree"] ]:
	s=subprocess.Popen( command, universal_newlines=True, stderr=subprocess.PIPE  )
	o,e = s.communicate()
	if s.returncode != 0:
		cleanup()
		sys.exit( f"Received error \"{e.strip()}\" running {command[0].split('/')[-1]}" )


#validate rearrangements output
#at the moment, it seems like IgPhyML output varies a bit even with a specified seed, so skip that test for now.
import hashlib
checksums = { "output/tables/tests_rearrangements.tsv":"896086eeac85cff23fd4c5924069a8fc" }#, "output/sequences/nucleotide/tests_inferredAncestors.fa":"46b03fc95afec5fda950f048d99db34b"}
for toValidate in checksums:
	h = hashlib.md5()
	with open(toValidate, 'rb') as handle:
		buf=handle.read()
		h.update(buf)
	if not h.hexdigest() == checksums[toValidate]:
		cleanup()
		sys.exit( f"{toValidate} failed validation!")


#All done!
cleanup()
print("All tests passed!")
