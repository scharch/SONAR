#!/usr/bin/env python3

"""
filterAIRR.py

This is a utility script to filter an AIRR-formatted TSV based on a set of 
    user-specified rules.

Usage: filterAIRR.py [options] RULE...

Options:
    --input INAIRR.tsv    An AIRR-formatted rearrangements file to be
                              [default: output/tables/<project>_rearrangements.tsv]
    --output OUTAIRR.tsv  Where to save the filtered rearrangements. Use
                              `stdout` to print to screen. [default: STDOUT]
    --or                  Flag to indicate that RULEs should be applied using
                              OR logic instead of AND. [default: False]
    RULE                  Filtering rules to apply. A RULE consists of three
                              parts, separated by commas:
                          1. COLUMN - the name of the column on which the rule
                              is to be applied.
                          2. OPERATOR - one of =~, !~, eq, ne, is, not, ==, !=, >=, 
                              >, <=, <, -, or !-. The first four assume that VALUE 
                              is a string; the next two treat COLUMN as a boolean;  
                              the next six will attempt to cast both COLUMN and 
                              VALUE as floats; and the last two  will treat VALUE 
                              as a file name from which to read a list of (eg) 
                              sequence_ids or cell_ids to be kept ('-') or 
                              discarded ('!-').
                          3. VALUE - the value to which COLUMN should be compared.
                              When enclosed in backticks, VALUE will be interpreted
                              as a second column name, allowing RULEs like
                              "n1_length,>,`n2_length`" to find rearrangements
                              with n1 insertions longer than their n2 insertions.
                              If OPERATOR is '-' or '!-' and VALUE is 'STDIN', then 
                              the list will be obtained from the streaming input.

Created by Chaim A Schramm on 2020-07-02.

Copyright (c) 2020 Vaccine Research Center, National Institutes of
                         Health, USA. All rights reserved.

"""

import airr, re, fileinput
import sys, os
from docopt import docopt

try:
	from SONAR import *
except ImportError:
	find_SONAR = sys.argv[0].split("SONAR/utilities")
	sys.path.append(find_SONAR[0])
	from SONAR import *


def processRule(sentence, reader):

	parts = sentence.split(",")

	if not len(parts) == 3:
		sys.exit(f"Invalid rule '{sentence}'\nRules must have three parts, separated by commas.")

	if not parts[0] in reader.fields:
		sys.exit(f"Error: Cannot find field {parts[0]} in {arguments['--input']}")

	valFormatted = f"\"{parts[2]}\""
	valTest = re.match("`(.+)`$",parts[2])
	if valTest:
		if not valTest.groups(0) in reader.fields:
			sys.exit(f"Error: Cannot find field {valTest.groups(0)} in {arguments['--input']}")
		else:
			valFormatted = f"r['{valTest.groups(0)}']"

	if parts[1] == "=~":
		rule = f"re.search({valFormatted}, r['{parts[0]}'])"

	elif parts[1] == "!~":
		rule = f"not re.search({valFormatted}, r['{parts[0]}'])"

	elif parts[1] == "eq":
		rule = f"r['{parts[0]}'] == {valFormatted}"

	elif parts[1] == "ne":
		rule = f"r['{parts[0]}'] != {valFormatted}"

	elif parts[1] == "is":
		rule = f"r['{parts[0]}']"

	elif parts[1] == "not":
		rule = f"not r['{parts[0]}']"

	elif parts[1] in ["==", "!=", ">=", ">", "<=", "<"]:
		rule = f"float(r['{parts[0]}']) {parts[1]} float({valFormatted})"

	elif parts[1] in ["-", "!-"]:
		if parts[2] == "STDIN":
			parts[2] = "-"
		toKeep = [ line.strip() for line in fileinput.input(parts[2]) ]
		rule = f"r['{parts[0]}'] in {toKeep}"
		if parts[1] == "!-":
			rule = f"r['{parts[0]}'] not in {toKeep}"

	else:
		sys.exit(f"Invalid operator {parts[1]}\nAllowed choices are =~, !~, eq, ne, is, not, ==, !=, >=, >, <=, <, -, or !-.")

	return rule


def main():

	reader = airr.read_rearrangement( arguments['--input'] )
	rules = [ processRule(rule, reader) for rule in arguments['RULE'] ]

	if arguments['--output'] != "STDOUT":
		sys.stdout = open(arguments['--output'], "w")

	writer = airr.io.RearrangementWriter( sys.stdout, fields=reader.fields )

	for r in filterAirrTsv(reader, rules, useOR=arguments['--or']):
		writer.write( r )



if __name__ == '__main__':

	arguments = docopt(__doc__)

	prj_tree = ProjectFolders(os.getcwd())
	prj_name = fullpath2last_folder(prj_tree.home)

	arguments['--input'] = re.sub("<project>", prj_name, arguments['--input'])
	if not os.path.isfile(arguments['--input']):
		sys.exit(f"Cannot find rearrangements file {arguments['--input']}")
	elif not airr.validate_rearrangement(arguments['--input']):
		sys.exit(f"File {arguments['--input']} is not in valid AIRR format.")

	#log command line
	logCmdLine(sys.argv)

	main()
