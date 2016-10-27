#!/usr/bin/env python

"""
splitFastqForUSearch.py

This script quickly splits large FASTQ files into sizes 
      that 32-bit USearch can handle.

Usage: splitFastqForUSearch.py <file.fastq> [ <read2.fastq> -o split -n 2000000 ]

Options: 
   -h --help             Show this documentation
   <file.fastq>          File to split
   <read2.fastq>         Optional second file for a coordinated split of 
                            un-merged paired-end reads
   -o --output SPLIT     Header for output file names, to which r1/r2 (if necessary),
                            a serial number, and .fq will be appended.
                            [default: split]
   -n --number 2000000   Number of reads in each output file [default: 2000000]

Created by Chaim A Schramm 2016-03-15
Copyright (c) 2016 Columbia University Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

"""

import sys
from docopt import docopt
from Bio.SeqIO.QualityIO import FastqGeneralIterator

try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *


def main():

    count = 0
    fnum  = 1

    if not paired:

        handle = open( "%s_%03d.fq" % (arguments['--output'], fnum), "w")
        for t,s,q in FastqGeneralIterator(open(arguments['<file.fastq>'], "r")):
            if count >= arguments['--number']:
                handle.close()
                count = 0
                fnum += 1
                handle = open( "%s_%03d.fq" % (arguments['--output'], fnum), "w")
            handle.write( "@%s\n%s\n+\n%s\n" % (t,s,q) )
            count += 1
                
    else:
        
        #going to assume all reads are in both files and skip error checking
        h1 = open( "%s_R1_%03d.fq" % (arguments['--output'], fnum), "w")
        h2 = open( "%s_R2_%03d.fq" % (arguments['--output'], fnum), "w")
        r2_gen = FastqGeneralIterator(open(arguments['<read2.fastq>'], "r"))

        for t,s,q in FastqGeneralIterator(open(arguments['<file.fastq>'], "r")):
            if count >= arguments['--number']:
                h1.close()
                h2.close()
                count = 0
                fnum += 1
                h1 = open( "%s_R1_%03d.fq" % (arguments['--output'], fnum), "w")
                h2 = open( "%s_R2_%03d.fq" % (arguments['--output'], fnum), "w")
            h1.write( "@%s\n%s\n+\n%s\n" % (t,s,q) )
            h2.write( "@%s\n%s\n+\n%s\n" % r2_gen.next() )
            count += 1



if __name__ == '__main__':

    saveSysArg = []+sys.argv

    arguments = docopt(__doc__)
    arguments['--number'] = int(arguments['--number'])

    #log command line
    logCmdLine(saveSysArg)

    paired = False
    if arguments['<read2.fastq>'] is not None:
        paired = True

    main()
