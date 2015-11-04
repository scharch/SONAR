#!/usr/bin/env python

"""
checkClusterBlast.py

This script monitors the SGE array jobs submitted by scripts 1.1 and 1.2 of the
      SONAR pipeline. When finished, it checks for failed BLAST jobs and
      resubmits as necessary. (Failure is defined as absence of a BLAST output
      file; there may be other failure modes that will not be caught.) If all 
      jobs completed successfully, can optionally start the next step of the
      pipeline.

Intended for internal use within SONAR only...

Usage: checkClusterBlast.py -gene <v|j|d|c> -big 100 -check check.sh
                            [ -after "next_step.py -args ..." ]

    Invoke with -h or --help to print this documentation.

    gene        Which gene blast should be monitored?
    big         Maximum index of the split files
    check       QSub script submitting this script to the cluster to continue
                    to monitor further rounds.
    after	Optional script (with all arguments) to be called after all
                    BLAST jobs have completed successfully.

Created by Chaim A Schramm on 2015-07-30.

Copyright (c) 2011-2015 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
find_SONAR_on_cluster = sys.argv[0].split("sonar/utilities")
sys.path.append(find_SONAR_on_cluster[0])
from sonar import *



def main():

    global gene, maxInd, myName, nextCmd

    getFolder   = dict( v="v",j="j",c="j",d="j" ) #since c and d are done in the same working folder as j
    addIn       = dict( v="",j="",c="_C",d="_D" ) #since c and d have slightly different file names

    folder      = "%s/%sgene"     % ( prj_tree.annotate, getFolder[gene] )
    mainProgram = "%s/%sblast.sh" % ( folder, gene )

    bad = 0
    for ind in range(maxInd):
        if os.path.isfile("%s/%s%s_%03d.txt"%(folder,prj_name,addIn[gene],ind+1)):
            pass
        else:
            subprocess.call(["qsub", "-t", str(ind+1), mainProgram])
            bad += 1

    if bad == 0:
        if nextCmd != "":
            subprocess.call(nextCmd.split(" "))
        print "All %s blast jobs finished successfully" % gene #goes to log file when called by monitor script
    else:
        subprocess.call(["qsub", myName])



if __name__ == '__main__':

	#check if I should print documentation
	q = lambda x: x in sys.argv
	if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]) or len(sys.argv)<7:
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, gene="gene", big="maxInd", check="myName", after="nextCmd")
	gene, maxInd, myName, nextCmd = getParasWithDefaults(dict_args, dict(nextCmd=""), "gene", "maxInd", "myName", "nextCmd")

        if not gene in ["v","j","c","d"]:
            print "gene %s not recognized. options are v d j c" % gene
            sys.exit(1)
            
        if not os.path.isfile(myName):
            print "Can't find self resubmission script %s" % myName
            sys.exit(1)

	prj_tree  = ProjectFolders(os.getcwd())
	prj_name  = fullpath2last_folder(prj_tree.home)

	main()
