#!/usr/bin/env python3

"""
This is the new SONAR setup/install script. It will check for
    prerequisites and configure some helper files.

Usage: setup.py
"""

import os,subprocess,sys,glob

SONAR_HOME=os.getcwd()
if len(glob.glob("%s/commonVars.py"%SONAR_HOME))==0:
    SONAR_HOME = sys.argv[0]
    if not os.path.isabs(SONAR_HOME):
        sys.exit("Can't find full path to SONAR home directory. You may need to call setup.py from within the SONAR directory or use the full absolute path.")

if not sys.platform.startswith("linux") and not sys.platform.startswith("darwin"):
    sys.exit("Error, cannot recognize OS. Expected 'linux' or 'darwin' (macos), but got '%s'"%sys.platform)

try:
    from docopt import docopt
except ImportError:
    sys.exit("DocOpt is a required library for SONAR. Please run `pip3 install docopt --user`")

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython is required for SONAR. Please run `pip3 install Biopython --user`")

try:
    import airr
except ImportError:
    sys.exit("AIRR is a required library for SONAR. Please run `pip3 install airr --user`")

try:
    from fuzzywuzzy import fuzz
except ImportError:
    print("fuzzywuzzy is not installed - the master script will not work.\nYou can fix this later by running `pip3 install fuzzywuzzy --user`.\nProceeding with install...\n\n",file=sys.stderr)

try:
    from ete3 import *
except ImportError:
    print("ete3 is not installed - tree plotting will not work.\nYou can fix this later by running `pip3 install ete3 --user`.\nProceeding with install...\n\n",file=sys.stderr)

try:
    from PyQt4.QtGui import QGraphicsSimpleTextItem, QGraphicsEllipseItem, QColor, QFont, QBrush, QPen
except ImportError:
    print("PyQt4 is not installed - tree plotting will not work.\nYou can fix this later by running `sudo apt-get install python-numpy python-qt4 python-lxml python-six`.\nProceeding with install...\n\n",file=sys.stderr)

try:
    import pandas
except ImportError:
    print("pandas is not installed - comparison of GSSPs (5.4) will not work.\nYou can fix this later by running `pip3 install pandas --user`.\nProceeding with install...\n\n",file=sys.stderr)

check = subprocess.call(["perl", "-MBio::SeqIO", '-e', '1'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if check == 1:
    sys.exit("BioPerl is a required for SONAR. Please run `cpanm Bio::Perl`")

check = subprocess.call(["perl", "-MList::Util", '-e', '1'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if check == 1:
    sys.exit("List::Util is a required library for SONAR. Please run `cpanm List::Util`")

check = subprocess.call(["perl", "-MAlgorithm::Combinatorics", '-e', '1'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if check == 1:
    sys.exit("Algorithm::Combinatorics is a required library for SONAR. Please run `cpanm Algorithm::Combinatorics`")

check = subprocess.call(["perl", "-MPDL::LinearAlgebra::Trans", '-e', '1'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if check == 1:
    sys.warn("PDL::LinearAlgebra::Trans is not installed - ancestor inference will not work.\nYou can fix this later by running `cpanm PDL::LinearAlgebra::Trans`.\nProceeding with install...\n\n")


#R library checks
for lib in ["docopt","ggplot2","MASS","grid"]:
    s=subprocess.Popen(['R','--vanilla','--slave','-e', '"%s" %%in%% installed.packages()[,"Package"]'%lib],
                       stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)
    o,e = s.communicate()
    if o.strip().split(" ")[1] == "FALSE":
        sys.exit("R Package %s is not installed. Please start R and run the command `install.packages('%s')`"%(lib,lib))

#cluster?
cluster_exists = ""
while cluster_exists.upper() not in ["Y", "N"]:
    cluster_exists = input("Is there a cluster available to use with SONAR [y/N]? ")
    if cluster_exists == "":
        cluster_exists = "N"

if cluster_exists.upper() == "Y":
    qsub = input("Please enter the command used to submit jobs to the cluster [qsub]: ")
    if qsub == "":
        qsub = "qsub"

#print out sonar, paths.py, and PPvars.pm
##################################################################
with open("%s/PPvars.pm"%SONAR_HOME, "w") as ppvars:
    ppvars.write("""#!/usr/bin/env perl

package PPvars;
use strict;

use vars '@ISA', '@EXPORT', '$NAME', '$VERSION', '$DATE', '$AUTHOR';
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(ppath);


sub ppath{
    return '%s/third-party/';
}

1;
""" % SONAR_HOME)
##################################################################
print_cluster = "clusterExists = False"
if cluster_exists.upper() == "Y":
    print_cluster = "clusterExists = True\nqsub = '%s'" % qsub

blast    = "blastn_linux64"
clustalo = "clustalo"
clustalw = "clustalw2"
muscle   = "muscle"
vsearch  = "vsearch"

if sys.platform.startswith("darwin"):
    blast    = "blastn_macos"
    clustalo = "clustalo_macos"
    clustalw = "clustalw2_macos"
    muscle   = "muscle_macos"
    vsearch  = "vsearch_macos"
##################################################################
with open("%s/paths.py"%SONAR_HOME, "w") as paths:
    paths.write("""
SCRIPT_FOLDER  = '%s'
blast_cmd = '%s/third-party/%s'
clustalo  = '%s/third-party/%s'
clustalw  = '%s/third-party/%s'
muscle    = '%s/third-party/%s'
vsearch   = '%s/third-party/%s'
%s
""" % (SONAR_HOME, SONAR_HOME, blast, SONAR_HOME, clustal, SONAR_HOME, clustalw, SONAR_HOME, muscle, SONAR_HOME, vsearch, print_cluster))
##################################################################
with open("%s/sonar"%SONAR_HOME, "w") as sonar:
    sonar.write("""#!/usr/bin/env python3

\"\"\"
sonar

This is a master script to allow easy access to SONAR scripts without
    needing to remember the exact commands or to add multiple
    directories to the path.

Usage: sonar COMMAND [ARGS...]

Options:
   COMMAND  Name of the SONAR script to run. Partial matches will be
                honored. In case of ambiguity, the program will exit
                with a list of possible matches.
   ARGS     Arguments to be passed to the script.

Created by Chaim A Schramm 2018-11-15.

Copyright (c) 2018 Vaccine Research Center, National Institutes
                      of Health, USA. All rights reserved.

\"\"\"

from docopt import docopt
import glob, os, subprocess, sys
from fuzzywuzzy import fuzz,process

def main():
	script_list  = [fn for fn in glob.glob(\"%s/*/*.py\") if not os.path.basename(fn).startswith(\"_\")] + glob.glob(\"%s/*/*.pl\") + glob.glob(\"%s/*/*.R\")

	match_script = process.extract(arguments['COMMAND'], script_list, limit=5, scorer=fuzz.partial_ratio)

	if match_script[0][1] == 100 and match_script[1][1] < 100:
		subprocess.call( [ match_script[0][0] ] + arguments['ARGS'] )
	else:
		print(\"Input program '%%s' is unclear. Did you mean one of the following?\"%%arguments['COMMAND'])
		for i in match_script:
			print(\"\\t\"+os.path.basename(i[0]))
		sys.exit()

if __name__ == '__main__':

	arguments = docopt(__doc__, options_first=True, version=\"SONAR v3.0\")
	main()

""" %(SONAR_HOME, SONAR_HOME, SONAR_HOME) )
##################################################################
os.chmod("%s/sonar"%SONAR_HOME, 0o755)
