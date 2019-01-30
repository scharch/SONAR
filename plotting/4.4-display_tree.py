#!/usr/bin/env python3

"""
4.4-plot_tree.py

This script uses the ete module to display figure-quality phylogenetic "birthday" trees from longitudinal data.
      Drawing engine requires an X server to run, use xvfb-run if needed.

Usage: 4.4-plot_tree.py -t newick.txt -n natives.csv
                        [ -o tree.png -c collapse.txt -i intermediates.csv -a annotations.txt --dates formattedTimes.txt ] 
                        [ --colors #FF000 ]... [ --sc 5000 | -w 4 ] [ --sp 12 | --he 4 ] [ -f 1 -r 300 ] 
                        [ --left --showAll --path --noDots --noUCA --noV --noGuide ]

Options:
    -t newick.txt                Tree to display, in Newick format. Leaf sequences are expected to fall into one of three
                                     categories: Germline/outgroup, whose label begins with "IG", "VH", "VK", or "VL";
                                     native antibodies, as identified in a tab-delimited text file provided with the -n
                                     option; or time-point-labeled sequences whose labels begins with on of the time point
                                     codes specified in the tab-delimited text file provided with the -n option. 
                                     Internal/non-terminal nodes can be labeled in any fashion (including not at all) and
                                     are ignored unless the -i or -c options are used (see below).
    -n natives.csv               Tab-delimited text file with the time points present in the tree (first column, in temporal 
                                     order), the sequences to be labeled as natives (second column), and optionally, 
                                     formatted labels to be displayed for them (third column). **ALL** time points must be 
                                     included in this file, even if no native antibodies have been identified from a 
                                     particular time point. In addition, time point labels must be fully unique. (No label 
                                     can be a prefix of another label; use 012 and 122 instead of 12 and 122. Of course, 
                                     this must still match what is in the newick file.) An example of the file format is 
                                     given below.
    -o tree.png                  File name for saving output image. PNG is recommended, but the ete rendering engine will
                                     recognize a .pdf extension (and possibly others) and respond accordingly. 
                                     [default: tree.png]
    -c collapse.txt              Optional tab-delimited text file specifying a list of internal nodes which should be 
                                     collapsed in the final display. First column is the id of the node to collapse (must be
                                     labeled in the newick file), second column is an optional message/label to print at the
                                     collapsed node. If no label is provided, the collapsed node will be annotated with the
                                     number of leaves hidden.
    -i intermediates.csv         Optional tab-delimited text file specifying intermediate nodes to highlight. First column
                                     is the node label in the newick file, optional second column gives a formatted label to
                                     display on the tree.
    -a annotations.txt           Optional tab-delimited text file specifying annotations for some or all nodes. First column
                                     is node label in the newick file, second is the annotation (eg CDR3 sequence). If the 
                                     annotation is "shapeXX", adds a shape annotation. First digit indicates circle (0) or 
                                     square (1); second digit is color (white, black, red, blue, green).
    --dates formattedTimes.txt   Optional tab-delimited text file to convert timepoint labels as used in the tree file and 
                                     the natives file (first column) to something more intelligible for the figure legend 
                                     (second column).
    --colors #FF000              Manually specify RGB hex colors for each time point. If not used, the program defaults to
                                     choosing up to 30 equally-ish spaced colors, sorted by hue from red (early) to blue (late).
    --sc 5000                    Display scale for the tree. According to the ete documentation, this value corresponds 
                                     to the size in pixels of a tree branch with unit length. Mutually exclusive with -w. 
                                     [default: 5000]
    -w 4                         Width (in inches) of the output image. Includes a quarter-inch margin around all borders and 
                                     a half-inch strip on the right side dedicated to a color guide, so the minimum value is 2 
                                     (1.5 if --noGuide is used). The program uses this value as a guideline to estimate a value
                                     for --sc, so the result will not be exact.
    --sp 12                      Vertical spacing of adjacent tree brances, in pixels. [default: 12]
    --he 4                       Height (in inches) of the output image. Includes a quarter-inch margin around all borders. 
                                     The program uses this value as a guideline to estimate a value for --sp, so the result
                                     will not be exact. Because spacing cannot be less than 0 and the program tries to ensure 
                                     that branches and labels are legible, there is a minimum possible height, and the program 
                                     will throw an error if the value of -h is less than that minimum. In the event that a 
                                     smaller tree is required, I recommend increasing -f and then reducing the output image
                                     as necessary. This will at at least keep the labels legible even as the branches start 
                                     to run togethere
    -f 1                         Magnification factor for text labels. Default size is approximately 0.1 inches high at the 
                                     specified resolution. (Font face is hard-coded as Arial.) [default: 1]
    -r 300                       Resolution (in dpi) of the output image. Please see note for --sc. [default: 300]
    --left                       Flag indicating that tree should be displayed facing left. [default: False]
    --showAll                    Flag indicating that labels should be displayed for all leaves. [default: False]
    --path                       Flag indicating that the evolutionary pathway(s) from root to the native antibody/ies should 
                                     be highlighted with a thicker line (does not change any colors). [default: False]
    --noDots                     Flag indicating that dots should not be used to highlight the postion(s) of native antibody/ies.
                                     [default: False]
    --noUCA                      Flag indicating that the root of the tree should not be labeled as "UCA". [default: False]
    --noV                        Flag indicating that the germline V sequence should not be labeled. [default: False]
    --noGuide                    Flag indicating that color guide should be hidden. [default: False]



Sample format of natives.csv:

039                                   #timepoints without natives (in order)
049        
060     VRC26m      CAP256-VRC26.1    #natives with nicely formatted display labels
120     VRC26       CAP256-VRC26.8
175        
206     VRC26k                        #not renamed, displayed as is


Created by Chaim A Schramm 2013-06-18
Edited and commented for publication by Chaim A Schramm on 2015-11-04 --happy birthday, Lisa <3
Edited to add showAll and noGuide options 2016-08-11 by CAS
Added left-facing option 2017-03-21 by CAS
Modified to add annotation option 2017-07-11 by CAS
Added date option 2017-09-27 by CAS
Modified to adjust automatic color selection and allow annotation of internal nodes 2017-09-28 by CAS
Modified to add shapes as annotations to branches 2018-04-12 by CAS
Edited to use Py3 and DocOpt by CAS 2018-08-29.

Copyright (c) 2013-2018 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

"""

import sys, os, re, colorsys
from docopt import docopt
from ete3 import *
from PyQt4.QtGui import QGraphicsSimpleTextItem, QGraphicsEllipseItem, QColor, QFont, QBrush, QPen
import string, random

try:
	from sonar.plotting import *
except ImportError:
	find_SONAR = sys.argv[0].split("sonar/plotting")
	sys.path.append(find_SONAR[0])
	from sonar.plotting import *


def parseNatives(infile):
	
	nats = dict()
	times = []

	reader	= csv.reader(open(infile, "r"), delimiter = sep)
	for row in reader:
		if len(row )== 3:
			date, name, show = row
		elif len(row) == 2:
			date, name = row
			show = name
		elif len(row) == 1:
			[ date ] = row
			name = ''
		else:
			print( "Unrecognized row while parsing natives from %s: %s" % (infile, row) )

		if date not in times:
			times.append(date)
		if name != '':
			nats[name] = dict( display=show, timepoint=date )

	if len(times) > len(arguments['--colors']):
		sys.exit( "Only %d colors specified. Please use the -colors options to input enough color selections for each of the %d time points detected.\n"%(len(arguments['--colors']),len(times)) )

	return nats, times


def parseInternalNodes(infile):
	if ( infile is None ):
		return []

	iNodes = dict()

	reader	= csv.reader(open(infile, "r"), delimiter = sep)
	for row in reader:
		if len(row)==2:
			name, show = row
		elif len(row)==1:
			[ name ] = row
			show	 = name
		else:
			print( "Unrecognized row while parsing node labels from %s: %s" % (infile, row) )

		if show == '':
			show = name

		iNodes[name] = show

	return iNodes


def parsePrettyDates(infile, unformatted):
	if ( infile is None ):
		return dict( zip(unformatted, unformatted) )

	dates	= dict()

	reader	= csv.reader(open(infile, "r"), delimiter = sep)
	for row in reader:
		if len(row)==2:
			name, show = row
		elif len(row)==1:
			[ name ] = row
			show	 = name
		else:
			print( "Unrecognized row while parsing dates from %s: %s" % (infile, row) )

		if show == '':
			show = name

		dates[name] = show

	return dates


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
		ns['shape'] = 'diamond'
		
		label = collapseList[node.name]
		if label == node.name:
			label = "Collapsed %d leaves in %d levels" % (len(node.get_leaves()), node.get_farthest_leaf(topology_only=True)[1]+1)
		tf = TextFace(" %s"%label,fgcolor="#B0B0B0",ftype='Arial',fsize=fontSize*.75, fstyle='italic')
		faces.add_face_to_node(tf, node, 0, position='branch-right')
		if not arguments['--noDots']: ns['size'] = 2 * ns['hz_line_width']
		ns['fgcolor'] = '#B0B0B0'
	else:
		#default appearance
		ns['size'] = 0
		ns['fgcolor'] = myRGB
		ns['vt_line_color'] = myRGB
		ns['hz_line_color'] = myRGB
		ns['hz_line_width'] = arguments['-r']/72
		ns['vt_line_width'] = arguments['-r']/72
		
		#highlight pathways to natives, if desired
		if arguments['--path'] and node.onPathway:
			ns['hz_line_width'] = arguments['-r']/24
			ns['vt_line_width'] = arguments['-r']/24
			
		#label natives
		if node.isNat:
			tf = TextFace(" %s"%node.name,fgcolor=myRGB,ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-right')
			if not arguments['--noDots']: ns['size'] = 2 * ns['hz_line_width']
			ns['fgcolor'] = '#000000'

		#label germline outgroup, if desired
		elif re.match("(IG|VH|VK|VL)", node.name) and not arguments['--noV']:
			tf = TextFace(" %s"%node.name,ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-right')
			ns['fgcolor'] = '#000000'

		#label all leaves, if desired
		elif arguments['--showAll'] and node.is_leaf():
			tf = TextFace(" %s"%node.name,ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-right')

		#add annotation, if provided
		if "annotation" in node.features:
			ann = re.findall("shape([01])([0-4]?)", node.annotation)
			if ann:
				#first digit is shape: 0=circle, 1=square (native to ete; other shapes require messing with Qt)
				#second digit is fill: 0/missing=white, 1=black, 2=red, 3=blue, 4=green
				fillColors = ['white','black','red','blue','green']
				whichCol = 0
				for shape in ann:
					thisFill   = fillColors[ int(shape[1] or 0) ]
					whichCol += 1
					if shape[0] == '0':
						faces.add_face_to_node(TextFace("  "), node, whichCol)
						whichCol += 1
						cf = CircleFace(fontSize/2, thisFill)
						faces.add_face_to_node(cf, node, whichCol)
					else:
						faces.add_face_to_node(TextFace("  "), node, whichCol)
						whichCol += 1
						rf = RectFace(fontSize-1, fontSize-1, 'black', thisFill)
						faces.add_face_to_node(rf, node, whichCol)
			else:
				tf = TextFace(" %s"%node.annotation,ftype='Arial',fsize=fontSize)
				if node.is_leaf():
					faces.add_face_to_node(tf, node, 0, position='aligned')
				else:
					faces.add_face_to_node(tf, node, 0, position='branch-bottom')				     

		#label intermediates, if desired
		if node.isIntermediate:
			f=faces.DynamicItemFace(iLabel, node.name)
			faces.add_face_to_node(f, node, 0, position='branch-right')

		#label UCA/root, if desired
		if not arguments['--noUCA'] and not node.up: 
			tf = TextFace(' UCA ',ftype='Arial',fsize=fontSize)
			faces.add_face_to_node(tf, node, 0, position='branch-top')
			# change horizontal line color??

	node.set_style(ns)


def makeRainbow(thisLevel, numLevels):
	if numLevels==1:
		return "#000000" # might need to check for a manual color input here
	else:
		subset = [int( a * (len(arguments['--colors'])-1) / (numLevels-1) ) for a in range(numLevels)]
	return arguments['--colors'][subset[thisLevel]]


def iLabel(node, *args, **kargs):

	#code for making specialized faces for intermediates mostly cribbed from the ete2 website example (though not interactive):
	# http://pythonhosted.org/ete2/tutorial/tutorial_drawing.html#creating-your-custom-interactive-item-faces

	my_label = node.name

	ellipse = QGraphicsEllipseItem(0,0,fontSize*2,fontSize*2) #I think the first two are coords of center; second pair is major/minor axis
	ellipse.setPen(QPen(QColor( 'black' )))
	ellipse.setBrush(QBrush(QColor( 'white' )))

	text = QGraphicsSimpleTextItem(my_label)
	text.setParentItem(ellipse)
	text.setBrush(QBrush(QColor("black")))
	font = QFont("Arial",fontSize*.9,weight=80)
	font.setLetterSpacing(1, 2) #add 2 pixels between letters for legibility
	text.setFont(font)

	#Center text according to masterItem size
	tw = text.boundingRect().width()
	th = text.boundingRect().height()
	center = ellipse.boundingRect().center()
	text.setPos(center.x()+1-tw/2, center.y()-th/2) #since the last letter has an extra 2 pixels after it from the spacing command, adjust center to compensate
    
	return ellipse






def main():

	global fontSize
	
	#read in tree
	myTree= Tree(arguments['-t'], format=1)
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
	intNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in intermediates), myTree.traverse() )
	for i in intNodes:
		if (not i.onPathway):
			print( "Warning: intermediate %s (%s) doesn't seem to be on a pathway to one of the natives...\n" % (i.name, intermediates[i.name]) )
		i.isIntermediate = True
		i.name = intermediates[i.name]

	#should we collapse any nodes?
	numCollapsed = 0
	collapseNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in collapseList), myTree.traverse() )
	for c in collapseNodes:
		c.collapse = True
		numCollapsed += len(c.get_leaves()) - 1 #will be wrong if one collapsed node is in the
							#subtree of another collapsed node, but that 
							#should be sufficiently rare that I am ignoring it

	#annotations?
	if arguments['-a'] is not None:
		reader	= csv.reader(open(arguments['-a'], "r"), delimiter = sep)
		for row in reader:
			if len(row)==2:
				try:
					(myTree & row[0]).add_feature( "annotation", row[1] )
				except coretype.tree.TreeError:
					print( "Can't annotate ", row[0], ": name not found!" )
			else:
				print( "Unrecognized annotation row: ", row )

				
	#for everything else, read off its timepoint
	for leaf in myTree.iter_leaves():
		if leaf.isNat or re.match("(IG|VH|VK|VL)", leaf.name):
			# these nodes are not expected to fit the standard time-labeling scheme
			continue
		leafTime = re.match(timeRegex,leaf.name)
		if leafTime is None:
			print( "Warning: Could not match leaf %s to one of the known time points (%s)" % (leaf.name,", ".join(timepoints)) )
		else:
			leaf.timepoint= leafTime.group()
		
	#now go through the timepoints ***in reverse order***
	#	   and mark each ancestor (for coloring purposes)
	for date in reversed(timepoints):
		for node in myTree.search_nodes(timepoint=date):
			while node.up:
				node.timepoint = date
				node = node.up
			
	#optimize graphical parameters
	estimatedH = ( len(myTree.get_leaves()) - numCollapsed )*(arguments['-r']/72) + fontSize*1.4*len(natives)	  #bare minimum based on line thickness and label size
	if arguments['--he'] is not None:
		pixelH = (float(arguments['--he'])-0.5) * arguments['-r']
		if estimatedH <= pixelH:
			arguments['--sp'] = (pixelH - estimatedH) / (len(myTree.get_leaves())-numCollapsed)		#easy, just increase spacing to match
		else:
			sys.exit("Specified image height of %2.1f is too small for the input tree. Minimum is ~%2.1f.\n"%(float(arguments['--he']), estimatedH/arguments['-r'] + .5) )

	if arguments['-w'] is not None:
		pixelW = (float(arguments['-w']) - 2) * arguments['-r'] #.5 for margin, .75 for color key, .5 for native labels (may not be necessary), .25 for branch behind UCA
		if arguments['--noGuide']: pixelW = (float(arguments['-w']) - 1.25)
		arguments['--sc'] = pixelW / myTree.get_farthest_leaf()[1]
		

	#overload package functions for nicer graphics
	from ete3.treeview import qt4_render
	qt4_render.add_scale  = _custom_add_scale
	qt4_render.add_legend = _custom_add_legend

	#basically done; set TreeStyle, call layout, and render!
	ts			  = TreeStyle()
	ts.layout_fn		  = layout
	ts.show_scale		  = True
	ts.show_leaf_name	  = False
	ts.branch_vertical_margin = arguments['--sp']
	ts.scale		  = arguments['--sc']
	ts.margin_left		  = arguments['-r']/4
	ts.margin_right		  = arguments['-r']/4
	ts.margin_bottom	  = arguments['-r']/4
	ts.margin_top		  = arguments['-r']/4

	if arguments['--left']:
		ts.orientation = 1
		
	#make color guide
	if not arguments['--noGuide']:
		for time in range( len(timepoints) ):
			color = makeRainbow(time, len(timepoints))
			bar = RectFace(arguments['-r']/4,arguments['-r']/24,color,color)
			bar.margin_top	  = arguments['-r']/24
			bar.margin_bottom = arguments['-r']/24
			bar.margin_right  = arguments['-r']/24
			bar.margin_left	  = arguments['-r']/24
			ts.legend.add_face(bar, column=0)
			text = TextFace(datesInLegend[timepoints[time]],ftype="Arial",fsize=fontSize*0.75,fgcolor="#000000")
			text.hz_align	  = 0
			text.vt_align	  = 1
			text.margin_left  = arguments['-r']/24
			text.margin_right = arguments['-r']/24
			ts.legend.add_face(text, column=1)

	myTree.dist=0.05

	myTree.render(arguments['-o'], dpi=arguments['-r'], tree_style=ts)

####################################

#Custom modifications of functions from within ete package
def _custom_add_scale(img, mainRect, parent):
			       
	from PyQt4 import QtGui
	from ete3.treeview.qt4_render import _EmptyItem

	length	    = arguments['-f'] * arguments['-r'] / 4
	length	    = img._scale * numpy.ceil( 100 * float(length) / img._scale ) / 100 #this is my OCD...
	length_text = float(length) / img._scale
    
	height	    = length/3
	scaleItem   = _EmptyItem()
	customPen   = QtGui.QPen(QtGui.QColor("black"), arguments['-r']/72)
    
	line = QtGui.QGraphicsLineItem(scaleItem)
	line2 = QtGui.QGraphicsLineItem(scaleItem)
	line3 = QtGui.QGraphicsLineItem(scaleItem)
	line.setPen(customPen)
	line2.setPen(customPen)
	line3.setPen(customPen)
    
	if arguments['--left']:
		length = -length

	line.setLine(0, height/2, length, height/2)
	line2.setLine(0, 0, 0, height)
	line3.setLine(length, 0, length, height)
	scale_text = "%0.2f" % (length_text)
	scale = QtGui.QGraphicsSimpleTextItem(scale_text)
	scale.setFont(QtGui.QFont("Arial", fontSize*0.75))
	scale.setParentItem(scaleItem)
	scale.setPos(length/3, -height)
	if arguments['--left']:
		scale.setPos(2*length/3, -height)
    
	scaleItem.setParentItem(parent)
    
	if arguments['--left']:
		scaleItem.setPos(mainRect.bottomRight())
		if not arguments['--noGuide']:
			from ete3.treeview.qt4_render import _FaceGroupItem
			scaleItem.moveBy(-_FaceGroupItem(img.legend,None).get_size()[0]+arguments['-r']/24, -img.margin_bottom)
		else:
			scaleItem.moveBy(-img.margin_right, -img.margin_bottom)
	else:
		scaleItem.setPos(mainRect.bottomLeft())
		scaleItem.moveBy(img.margin_left, -img.margin_bottom)
	mainRect.adjust(0, 0, 0, height)


def _custom_add_legend(img, mainRect, parent):

	from ete3.treeview.qt4_render import _FaceGroupItem
	
	legend = _FaceGroupItem(img.legend, None)
	legend.setup_grid()
	legend.render()
	lg_w, lg_h = legend.get_size()
	mid_height = mainRect.height() / 2
	place_at = mid_height - lg_h/2
	legend.setParentItem(parent)
	legend.setPos(mainRect.width(), place_at)
	#legend.moveBy(-img.margin_right,0)
	mainRect.adjust(0, 0, lg_w+img.margin_right, 0)
	

############################################

if __name__ == '__main__':

	if os.environ.get('DISPLAY') is None:
		sys.exit("\nThe Qt4 drawing engine requires an active X server.\nPlease enable X forwarding or use xvfb-run.\n")


	arguments = docopt(__doc__)

	arguments['-f']	  = float(arguments['-f'])
	arguments['-r']	  = int(arguments['-r'])
	arguments['--sc'] = int(arguments['--sc'])
	arguments['--sp'] = int(arguments['--sp'])


	if arguments['-w'] is not None and (float(arguments['-w']) < 2 or (arguments['--noGuide'] and  float(arguments['-w']) < 1.5)):
		sys.exit("Minimum figure width is 2 inches with color guide or 1.5 without\n");
	

	#convert font magnification to actual size
	fontSize = arguments['-f'] * 20 * arguments['-r']/300


	#colors generated stochastically via http://tools.medialab.sciences-po.fr/iwanthue/ and sorted for increasing hue
	if len(arguments['--colors']) == 0:
		arguments['--colors'] = ["#BE4229", "#E74721", "#8C431C", "#CB6A27", "#E98C25", "#946E13", "#D0A620", "#8B8A22", "#A9B71D", "#555C10", "#839A21", "#82B531", "#4D8426", "#55C42A", "#34992A", "#2B6B2E", "#4FC456", "#33B66D", "#4296CB", "#5A8CE2", "#3E5988", "#656CE2", "#524EA0", "#8F83CC", "#A57CE4", "#8E46AD", "#C056EB", "#CA6BE4", "#7B4D87", "#D186D7"]

	#load natives
	natives, timepoints = parseNatives(arguments['-n'])
	timeRegex = re.compile( "(" + ")|(".join(timepoints) + ")" )

	#formatted dates?
	datesInLegend = parsePrettyDates( arguments['--dates'], timepoints )
	
	#parse intermediates and nodes to collapse
	# (function returns empty list if file is not specified)
	intermediates = parseInternalNodes(arguments['-i'])
	collapseList  = parseInternalNodes(arguments['-c'])


	#log command line
	logCmdLine(sys.argv)


	main()

