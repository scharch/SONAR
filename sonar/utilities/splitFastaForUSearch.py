#!/usr/bin/env python

"""
splitFastaForUSearch.py

This script quickly splits large FASTA files into sizes 
      that 32-bit USearch can handle.

Usage: splitFastaForUSearch.py <file.fasta> [ -o split -n 2000000 ]

Options: 
   -h --help             Show this documentation
   <file.fasta>          File to split
   -o --output SPLIT     Header for output file names, to which a serial
                            number and .fa will be appended. [default: split]
   -n --number 3000000   Number of reads in each output file [default: 3000000]

Created by Chaim A Schramm 2016-03-15
Duplicated for fastA version 2016-05-18
Copyright (c) 2016 Columbia University Vaccine Research Center, National
                          Institutes of Health, USA. All rights reserved.

"""

from docopt import docopt
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
	from sonar import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/utilities")
	sys.path.append(find_SONAR[0])
	from sonar import *

        

def main():

    count = 0
    fnum  = 1
    
    handle = open( "%s_%03d.fa" % (arguments['--output'], fnum), "w")
    for t,s in SimpleFastaParser(open(arguments['<file.fasta>'], "r")):
        if count >= arguments['--number']:
            handle.close()
            count = 0
            fnum += 1
            handle = open( "%s_%03d.fa" % (arguments['--output'], fnum), "w")
        handle.write( ">%s\n%s\n" % (t,s) )
        count += 1
                

if __name__ == '__main__':
    
    saveSysArg = []+sys.argv
    
    arguments = docopt(__doc__)
    arguments['--number'] = int(arguments['--number'])

    #log command line
    logCmdLine(saveSysArg)

    main()
