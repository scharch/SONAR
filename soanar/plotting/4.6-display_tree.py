#!/usr/bin/env python
# encoding: utf-8
"""
prettyTree.py

Created by Chaim A Schramm 2013-06-18
Copyright (c) 2013 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: prettyTree.py -t newick.txt -n natives.csv -o output.pdf [-i intermediates.txt] [-v <3|6>]

natives.csv should be in two column format with name and timepoint. 
if there are timepoints with no natives, they should be included as a row with empty first column.
optional third column with printname. obviously if no third column, original name will be displayed.
so eg:

VRC26   120   CAP256-VRC26.8
VRC26k  206                     #not renamed
VRC26m  060   CAP256-VRC26.1
        039                     #other timepoints
        049
        175

intermediates are not implemented yet, but will probably be a list of node names (so yes, might have to insert manually into Newick first) and possibly a second column for renaming (so that '271' can be displayed as 'I4')

v is 'version' -should be either figure *3* or *6*

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
		if len(row)==3:
			name, date, show = row
		elif len(row)==2:
			name, date = row
			show = name
		else:
			print "Unrecognized row: ", row

		if date not in times:
			times.append(date)
		if show == '':
			show = name
		if name != '':
			nats[name] = dict( display=show, timepoint=date )

	times = sorted(times,reverse=True)
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



def layout(node):
	ns = NodeStyle()
	ns['size'] = 0
	#ns['hz_line_width']=2
	#ns['vt_line_width']=2

	#if node.up is None:
	#	print "ROOT"
	#else:
	#	print node.name

	colors = ['#990099', '#0000FF', '#00C8FF', '#339900', '#FF9900', '#FF0000']
	myRGB = colors[int(node.timepoint)] #makeRainbow(node.timepoint, len(timepoints))

	ns['fgcolor'] = myRGB
	ns['vt_line_color'] = myRGB
	ns['hz_line_color'] = myRGB
	
	
	if version==6 and not doNotCollapse(node):
		ns['hz_line_type'] = 1
		ns['vt_line_color'] = '#FFFFFF'
		ns["draw_descendants"] = False
		node.set_style(ns)
		return
        

	#if node.onPathway:
	#	ns['hz_line_width'] = 2
	#	ns['vt_line_width'] = 3

	if node.isNat:
		if node.name == "(to CAP256-VRC26.02-12)":
			myRGB='#000000' #kludge
			ns['hz_line_type'] = 1
		tf = TextFace(" %s"%node.name,fgcolor=myRGB,ftype='Arial',fsize=8)
		#faces.add_face_to_node(tf, node, 0, position='branch-right')
		if version==3:
			ns['size'] = 7
			ns['fgcolor'] = '#000000'

	if re.match("IG", node.name):
		ns['hz_line_color']='#FFFFFF'
		tf = TextFace(" %s"%node.name,ftype='Arial',fsize=8)
		if version==3:
			#faces.add_face_to_node(tf, node, 0, position='branch-right')
			ns['size'] = 7
			ns['fgcolor'] = '#000000'

	if node.isIntermediate:
		f=faces.DynamicItemFace(iLabel, node.name)
		faces.add_face_to_node(f, node, 0, position='branch-right')

	if not node.up: #UCA!
		tf = TextFace('UCA ',ftype='Arial',fsize=8)
                #if version==3:
		#faces.add_face_to_node(tf, node, 0, position='branch-top')
		#ns['hz_line_color']='#FF0000'
		#else:
		#	ns['hz_line_color']='#FFFFFF'


	node.set_style(ns)

def makeRainbow(current, total):
	fraction = 1 - float(current+1)/total
	sat = 1
	if fraction > 2.0/9: sat = .5
	a = [255 * i for i in colorsys.hsv_to_rgb( fraction, 1, 1)]
	#print a, "#%02x%02x%02x"%tuple(a)
	return "#%02x%02x%02x"% tuple(a)


def doNotCollapse(node):
	
	#if not node.is_leaf():
	#	return False

	if node.onPathway or (node.up is None):
		#print 'keep %s'%node.name
		return True

	for a in node.iter_ancestors():
		#print node.name, a.name, intermediates.keys()
		if a.name in intermediates.values():
			#print "%s comes from %s"%(node.name,a.name)
			return True
			
	#print "collapse %s"%node.name
	return False


def main():
	
	myTree= Tree(treefile, format=1)
	for node in myTree.traverse():
		node.add_features(isNat=False, isIntermediate=False, onPathway=False, timepoint=-1)

	#first, find the mAbs
	natNodes = filter( lambda x: any( re.match("%s$"%n, x.name) for n in natives), myTree.iter_leaves() )

	#tell the mAbs that's what they are and mark the pathway down from UCA
	for mAb in natNodes:
		mAb.isNat = True
		mAb.timepoint = timepoints.index(natives[mAb.name]['timepoint'])
		mAb.name = natives[mAb.name]['display']
		level = mAb
		while level.up:
			level.onPathway = True
			level = level.up

	#mark the intermediates
	for iNode in intermediates:
		ii = myTree.search_nodes(name=iNode)
		if (len(ii) == 0):
			print "Warning: no intermediate %s"%iNode
			continue
		i = ii[0]
		if (not i.onPathway):
			print "Warning: intermediate %s (%s) doesn't seem to be on pathway...\n" % (i.name, intermediates[i.name])
		i.isIntermediate = True
		i.name = intermediates[i.name]

	#for everything else, read off its timepoint
	for leaf in myTree.iter_leaves():
		#doNotCollapse(leaf)
		if leaf.isNat or re.match("IG", leaf.name):
			continue #they don't fit the scheme in the next line
		leaf.timepoint = timepoints.index(leaf.name[0:3]) #probably need a more flexible/general way of doing this
		
	#now that we can group them easily, go through the timepoints in reverse order and mark each ancestor (for coloring purposes)
	for index in range(len(timepoints)):
		for node in myTree.search_nodes(timepoint=index): #filter( lambda x: x.timepoint == index, myTree.get_leaves() ):
			while node.up:
				node.timepoint = index
				node = node.up
				#node.name = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(6))
	myTree.timepoint = len(timepoints)-1
			
	#basically done; set TreeStyle, call layout, and render!
	ts = TreeStyle()
	ts.layout_fn = layout
	ts.show_scale = True
	ts.show_leaf_name = False
	if version==6:
		ts.branch_vertical_margin = 2
		ts.scale = 400
	else:
		ts.scale = 1500
	#ts.force_topology = True
	#ts.mode='c'

	#print "\nfinal processing..."
	#myTree.write( format_root_node=False, format=1, is_leaf_fn=doNotCollapse, outfile='foo.tree' )
	#collapsed = Tree( myTree.write( features=[], format_root_node=True, format=1, is_leaf_fn=doNotCollapse ), format=1 )
	myTree.dist=0.1
	myTree.render(outfile, dpi=600, tree_style=ts, h=10.5, units="in")



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



if __name__ == '__main__':
	# get parameters from input
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, t="treefile", n="nativefile", o="outfile", i="ifile", v="version")
	treefile, nativefile, outfile, ifile, version = getParasWithDefaults(dict_args, dict_args, "treefile", "nativefile", "outfile", "ifile", "version")

	natives, timepoints = parseNatives(nativefile)
	intermediates = parseIntermediates(ifile)

	if (version!=3 and version!=6):
		version =3

	main()

