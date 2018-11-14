#!/usr/bin/env python3

"""
checkClusterBlast.py

This script monitors the SGE array jobs submitted by scripts 1.1 and 1.2 of the
      SONAR pipeline. When finished, it checks for failed BLAST jobs and
      resubmits as necessary. (Failure is defined as absence of a BLAST output
      file; there may be other failure modes that will not be caught.) If all 
      jobs completed successfully, can optionally start the next step of the
      pipeline.

Intended for internal use within SONAR only...

Usage: checkClusterBlast.py --gene v --big 100 --check check.sh [ --after next.py --rehold otherJID ]

Options:
    --gene v            Which gene blast should be monitored? (v, d, j, or c)
    --big 100           Maximum index of the split files
    --check check.sh    QSub script submitting this script to the cluster to continue
                            to monitor further rounds.
    --after next.py     Optional script (with all arguments, enclosed in single quotes)
                            to be called after all BLAST jobs have completed successfully.
    --rehold otherJID   Optional pending jobname on which to renew a hold, if this 
                            current job needs to be resubmitted

Created by Chaim A Schramm on 2015-07-30.
Added rehold option CAS 2016-01-13.
Edited to use Py3 and DocOpt by CAS 2018-08-28.

Copyright (c) 2011-2018 Columbia University and Vaccine Research Center, National
                         Institutes of Health, USA. All rights reserved.

"""

import sys
import subprocess
from docopt import docopt
try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *




def main():

	getFolder   = dict( v="v",j="j",c="j",d="j" ) #since c and d are done in the same working folder as j
	addIn	    = dict( v="",j="",c="_C",d="_D" ) #since c and d have slightly different file names

	folder	    = "%s/%sgene"     % ( prj_tree.annotate, getFolder[arguments['--gene']] )
	mainProgram = "%s/%sblast.sh" % ( folder, arguments['--gene'] )

	bad = 0
	for ind in range(arguments['--big']):
		if os.path.isfile("%s/%s%s_%03d.txt"%(folder,prj_name,addIn[arguments['--gene']],ind+1)):
			pass
		else:
			subprocess.call(["qsub", "-t", str(ind+1), mainProgram])
			bad += 1
		    
	if bad == 0:
		if arguments['--after'] is not None:
			subprocess.call( re.split("\s+",arguments['--after']) )
		print( "All %s blast jobs finished successfully" % arguments['--gene'] ) #goes to log file when called by monitor script
	else:
		p = subprocess.Popen(["qsub", arguments['--check']], stdout=subprocess.PIPE)
		output, err = p.communicate()
		print( err )
		print( output )
		jobName = output.split('"')[1]
		if arguments['--rehold'] is not None:
			#This does not currently work on the C2B2 HPC due to permissions/administrative settings haven't
			#	  been able to find a work-around, though, so leaving it here in case it works for others.
			#Without this, if some c/d blast jobs fail, the jmonitor calling 1.3 won't wait for the resubmitted
			#	  jobs to finish, which could then cause an infinite-resubmission loop after 1.3 deletes
			#	  the fasta files needed to run c/d blast.
			#In practice, though, this whole que-monitor structure seems to be overkill/unnecesarry now that I've
			#	  updated a lot of other code.
			subprocess.call(["qalter", "-hold_jid", jobName, arguments['--rehold']])



if __name__ == '__main__':

	arguments=docopt(__doc__)
	arguments['--big'] = int( arguments['--big'] )

	if not arguments['--gene'] in ["v","j","c","d"]:
		print( "gene %s not recognized. options are v d j c" % arguments['--gene'] )
		sys.exit(1)
	    
	if not os.path.isfile(arguments['--check']):
		print( "Can't find self resubmission script %s" % arguments['--check'] )
		sys.exit(1)

	
	#log command line
	logCmdLine(sys.argv)

	prj_tree  = ProjectFolders(os.getcwd())
	prj_name  = fullpath2last_folder(prj_tree.home)

	main()
