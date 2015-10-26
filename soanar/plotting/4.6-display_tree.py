#!/usr/bin/env python

"""
4.6-plot_tree.py

This script uses the ete2 module to display figure-quality phylogenetic 
      "birthday" trees from longitudinal data.

Usage: 4.6-plot_tree.py -t newick.txt -n natives.csv
                        [ -o tree.png -c collapse.txt -i intermediates.csv -f 40
		          -m 10 -r 600 -s 3600 -w 6 -path -noDots -noUCA -noV    ]

    Invoke with -h or --help to print this documentation.

    t         -    Tree to display, in Newick format. Leaf sequences are expected
                         to fall into one of three categories: Germline/outgroup,
			 whose label begins with "IG", "VH", "VK", or "VL";
			 native antibodies, as identified in a tab-delimited text
			 file provided with the -n option; or time-point-labeled 
			 sequences whose labels begins with on of the time point
			 codes specified in the tab-delimited text file provided
			 with the -n option. Internal/non-terminal nodes can be
			 labeled in any fashion (including not at all) and are 
			 ignored unless the -i or -c options are used (see 
			 below).
    n         -    Tab-delimited text file with the time points present in the 
                         tree (first column, in temporal order), the sequences
			 to be labeled as natives (second column), and 
			 optionally, formatted labels to be displayed for them
			 (third column). **ALL** time points must be included
			 in this file, even if no native antibodies have been
			 identified from a particular time point. In addition,
			 time point labels must be fully unique. (No label can
			 be a prefix of another label; use 012 and 122 instead 
			 of 12 and 122. Of course, this must still match what is
			 in the newick file.) An example of the file format is 
			 given below.
    o         -    File name for saving output image. PNG is recommended, but
                         the ete2 rendering engine will recognize a .pdf
			 extension (and possibly others) and respond accordingly.
			 Default = tree.png
    c         -    Optional single-column text file specifying a list of internal
                         nodes which should be collapsed in the final display.
    i         -    Optional tab-delimited text file specifying intermediate nodes
                         to highlight. First column is the node label in the
			 newick file, optional second column gives a formatted
			 label to display on the tree.
    f         -    Desired font size for text labels. (Font face is hard-coded as
                         Arial.) Default = 40.
    m         -    Desired vertical spacing of tree brances. ete2 documentation
                         has this as the number of pixels between adjacent
			 branches; please see note for the -s parameter.
			 Default = 10.
    r         -    Desired resolution of the output image. Please see note for 
                         the -s parameter. Default = 600 dpi.
    s         -    Desired display scale for the tree. According to ete2
                         documentation, this value corresponds to the length in
			 pixels of a tree branch with unit length. In practice,
			 the effects of this parameter are inconsistent and 
			 interact in important ways with the -r, -s, and -w
			 parameters. The defaults for these 4 should work well
			 for most applications, but in the event that one needs
			 to be changed, the other 3 will likely need to be
			 adjusted via trial-and-error, as well. Default = 3600.
    w         -    Desired width (in inches) of the output image. Together with
                         options -r, -s, and -m, this will also fix image height,
			 which is therefore not provided as a free parameter.
			 Please see note for the -s parameter. Default = 6.
    path      -    Flag indicating that the evolutionary pathway(s) from root to
                         the native antibody/ies should be highlighted with a
			 thicker line (does not change any colors). Default = No.
    noDots    -    Flag indicating that dots should not be displayed to indicate
                         the postion(s) of native antibody/ies.
			 Default = dots are displayed.
    noUCA     -    Flag indicating that the root of the tree should not be
                         labeled as "UCA". Default = "UCA" label is shown.
    noV       -    Flag indicating that the germline V sequence should not be 
                         labeled. Default = germline V is labeled.



Sample format of natives.csv:

039                                   #timepoints without natives (in order)
049        
060     VRC26m      CAP256-VRC26.1    #natives with nicely formatted display labels
120     VRC26       CAP256-VRC26.8
175        
206     VRC26k                        #not renamed, displayed as is


Created by Chaim A Schramm 2013-06-18
Edited and commented for publication by Chaim A Schramm on 2015-10-26.

Copyright (c) 2013-2015 Columbia University and Vaccine Research Center, National
                               Institutes of Health, USA. All rights reserved.

"""


import sys, os, re, colorsys
from mytools import *
from ete2 import *
from PyQt4.QtGui import QGraphicsSimpleTextItem, QGraphicsEllipseItem, QColor, QFont, QBrush
import string, random


def parseNatives(infile):
	nats = dict()
	times = []

	reader 	= csv.reader(open(infile, "rU"), delimiter = sep)
	for row in reader:
		if len(row )== 3:
			date, name, show = row
		elif len(row) == 2:
			date, name = row
			show = name
		elif len(row) == 1:
			date = row
		else:
			print "Unrecognized row: ", row

		if date not in times:
			times.append(date)
		if name != '':
			nats[name] = dict( display=show, timepoint=date )

	if len(times) > 30:
		sys.exit( "Sorry, program is currently only configured to handle up to 30 time points.\n" )

	return nats, times


def parseIntermediates(infile):
	if ( infile is None ):
		return []

        iNodes = dict()

        reader  = csv.reader(open(infile, "rU"), delimiter = sep)
        for row in reader:
                if len(row)==2:
                        name, show = row
                elif len(row)==1:
                        name = row
                        show = name
                else:
                        print "Unrecognized row: ", row

                if show == '':
                        show = name

		iNodes[name] = show

	return iNodes


def parseCollapse(infile):
	if ( infile is None ):
		return []

	with open(infile) as f:
		return [line.strip() for line in f.readlines()]


def layout(node):

	# Get birthday color (defaults to black in case of error)
	myRGB = "#000000"
	if node.timepoint is not None:
		myRGB = makeRainbow(timepoints.index(node.timepoint), len(timepoints))

	#set up style object
	ns = NodeStyle()

	if node.collapse:
		ns['hz_line_type'] = 1
		ns['vt_line_color'] = '#FFFFFF'
		ns["draw_descendants"] = False
	else:
		#default appearance
		ns['size'] = 0
		ns['fgcolor'] = myRGB
		ns['vt_line_color'] = myRGB
		ns['hz_line_color'] = myRGB
		
		#highlight pathways to natives, if desired
		if pathway and node.onPathway:
			ns['hz_line_width'] = 2
			ns['vt_line_width'] = 3
			
		#label natives
		if node.isNat:
			tf = TextFace(" %s"%node.name,fgcolor=myRGB,ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-right')
			if dots: ns['size'] = 7
			ns['fgcolor'] = '#000000'

		#label germline outgroup, if desired
		if vgene and re.match("(IG|VH|VK|VL)", node.name):
			ns['hz_line_color']='#FFFFFF'
			tf = TextFace(" %s"%node.name,ftype='Arial',fsize=fontSize)
			ns['fgcolor'] = '#000000'

		#label intermediates, if desired
		if node.isIntermediate:
			f=faces.DynamicItemFace(iLabel, node.name)
			faces.add_face_to_node(f, node, 0, position='branch-right')

		#label UCA/root, if desired
		if uca and not node.up: 
			tf = TextFace('UCA ',ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-top')
			# change horizontal line color??

	node.set_style(ns)


def makeRainbow(thisLevel, numLevels):
	#colors generated stochastically via http://tools.medialab.sciences-po.fr/iwanthue/ and sorted for increasing hue
	fullList = ["#BE4229", "#E74721", "#8C431C", "#CB6A27", "#E98C25", "#946E13", "#D0A620", "#8B8A22", "#A9B71D", "#555C10", "#839A21", "#82B531", "#4D8426", "#55C42A", "#34992A", "#2B6B2E", "#4FC456", "#33B66D", "#4296CB", "#5A8CE2", "#3E5988", "#656CE2", "#524EA0", "#8F83CC", "#A57CE4", "#8E46AD", "#C056EB", "#CA6BE4", "#7B4D87", "#D186D7"]
	subset = [int( a * 30 / numLevels ) for a in range(numLevels)]
	return fullList[subset[thisLevel]]


def iLabel(node, *args, **kargs):

	#code for making specialized faces for intermediates mostly cribbed from the ete2 website example (though not interactive):
	# http://pythonhosted.org/ete2/tutorial/tutorial_drawing.html#creating-your-custom-interactive-item-faces

	my_label = args[0][0] #or maybe just node.name?

	ellipse = QGraphicsEllipseItem(0,0,16,16) #I think the first two are coords of center; second pair is major/minor axis
	ellipse.setBrush(QBrush(QColor( 'black' )))

	text = QGraphicsSimpleTextItem(my_label)
	text.setParentItem(ellipse)
	text.setBrush(QBrush(QColor("white")))
	text.setFont(QFont("Arial",8))

	#Center text according to masterItem size
	tw = text.boundingRect().width()
	th = text.boundingRect().height()
	center = ellipse.boundingRect().center()
	text.setPos(center.x()-tw/2, center.y()-th/2)
    
	return ellipse






def main():
	
	#read in tree
	myTree= Tree(treeFile, format=1)
	for node in myTree.traverse():
		node.add_features(isNat=False, isIntermediate=False, onPathway=False, collapse=False, timepoint=None)

	#first, find the mAbs
	natNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in natives), myTree.iter_leaves() )

	#tell the mAbs that's what they are and mark the pathway down from UCA
	for mAb in natNodes:
		mAb.isNat = True
		mAb.timepoint = natives[mAb.name]['timepoint']
		mAb.name = natives[mAb.name]['display']
		level = mAb
		while level.up:
			level.onPathway = True
			level = level.up

	#mark the intermediates
	intNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in intermediates), myTree.iter_leaves() )
	for i in iNodes:
		if (not i.onPathway):
			print "Warning: intermediate %s (%s) doesn't seem to be on a pathway to one of the natives...\n" % (i.name, intermediates[i.name])
		i.isIntermediate = True
		i.name = intermediates[i.name]

	#should we collapse any nodes?
	collapseNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in collapseList), myTree.iter_leaves() )
	for c in collapseNodes:
		c.collapse = True

	#for everything else, read off its timepoint
	for leaf in myTree.iter_leaves():
		if leaf.isNat or re.match("(IG|VH)", leaf.name):
			# these nodes are not expected to fit the standard time-labeling scheme
			continue 
		leaf.timepoint = re.match(timeRegex,leaf.name.group())
		if leaf.timepoint is None:
			print "Warning: Could not match leaf %s to a known time point!\n" % leaf.name
		
	#now go through the timepoints ***in reverse order***
	#	   and mark each ancestor (for coloring purposes)
	for date in reversed(timepoints):
		for node in myTree.search_nodes(timepoint=date):
			while node.up:
				node.timepoint = date
				node = node.up
			
	#basically done; set TreeStyle, call layout, and render!
	ts = TreeStyle()
	ts.layout_fn = layout
	ts.show_scale = True
	ts.show_leaf_name = False

	#if we set these, can only set either width OR height
	ts.branch_vertical_margin = margin
	ts.scale = scale

	myTree.dist=0.05
	myTree.render(outFile, dpi=res, tree_style=ts, w=width, units="in")




if __name__ == '__main__':
	# get parameters from input
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)


	#check various flags
	q = lambda x: x in sys.argv
	
	#hide dots at natives?
	dots = True
	if q("-noDots"):
		sys.argv.remove("-noDots")
		dots = False

	#bold pathways to natives
	pathway = False
	if q("-path"):
		sys.argv.remove("-path")
		pathway = True

	#don't label UCA/root?
	uca = True
	if q("-noUCA"):
		sys.argv.remove("-noUCA")
		uca = False

	#don't label Vgene/outgroup?
	vgene = True
	if q("-noV"):
		sys.argv.remove("-noV")
		vgene = False


	#parse remainder of options
	dict_args = processParas(sys.argv, t="treeFile", n="nativeFile", o="outFile", i="intFile",
				 c="collapseFile", s="scale", m="margin", r="res", f="fontSize", w="width")
	dict_defaults = dict( scale=3600, margin=10, res=600, fontSize=40, width=6, outFile="tree.png" )
	treeFile, nativeFile, outFile, intFile, collapseFile, scale, margin, res, fontSize, width = \
	    getParasWithDefaults( dict_args, "treeFile", "nativeFile", "outFile", "intFile", "collapseFile", "scale", "margin", "res", "fontSize", "width" )

	#load natives
	natives, timepoints = parseNatives(nativeFile)
	timeRegex = re.compile( "(" + ")|(".join(timepoints) + ")" )

	#parse intermediates and nodes to collapse
	# (functions return empty list if file is not specified)
	intermediates = parseIntermediates(intFile)
	collapseList = parseCollapse(collapseFile)


	main()

