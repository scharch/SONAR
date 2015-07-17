#!/usr/bin/env python
# encoding: utf-8
"""
1.3-heatmap.py

Created by Zhenhai Zhang on 2011-06-02.
Copyright (c) 2011 Columbia University Vaccine Research Center, National Institutes of Health, USA. All rights reserved.

Usage: 1.4-heatmap.py -l is_light [-r list_of_reads.txt] [-t title]


"""
import sys, os
from mytools import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def main():
	
	infiles = glob.glob("%s/*_divid.tab" %prj_tree.clustal_data)
	#print infiles
	
	for infile in infiles:
		print "processing file: %s...." %infile
		f_head = infile[ infile.rindex("/") + 1 : infile.rindex("_")]
		print "retrieving diverences and identities..."
		reader 	= csv.reader(open(infile, "rU"), delimiter = sep)
		natives = reader.next()[2:]
		
		# list contain dictionaries, each one corresponds to one native
		divids	= []
		for i in range(len(natives)):
			divids.append(dict())

		for row in reader:
			if len(row) - 2 < len(natives):
				pass
			elif len(r_list)>0 and (row[0] not in r_list):
				pass #simple, if not necessarily elegant
			else:
				read_id, divergence, identities = row[0], float(row[1]), map(float, row[2:])
				for ind, identity in enumerate(identities):
					divids[ind][read_id] = (divergence, identity)

		print "preparing array..."
		
		nat_divid = zip(natives, divids)
		
		for native, divid in nat_divid:
			print native
			
			my_array = generate_array(divid)
			my_array = manipulate_array(my_array)

			fig 	= plt.figure(figsize=(2, 2), dpi=600)
			
			ax	= fig.add_subplot(111)
			
			# adjust the figure to fill out the space.
			#adjustprops = dict(left=0.05, bottom=0.005, right=1.0, top=0.95)
			#fig.subplots_adjust(**adjustprops)
			
			#ax.set_title(fig_title, fontsize=8, fontweight="bold")

			#plt.contourf(my_array)
			ax.contour(my_array, origin='lower')
			ax.set_title("2D Plot of %s Reads"%title, fontsize=10)

			plt.xlim(xmin = 0, xmax = 40)
			plt.ylim(ymin = 60, ymax = 100)
			ax.get_yaxis().set_ticks([60,70,80,90,100])
			ax.get_xaxis().set_ticks([0,10,20,30,40])
			plt.setp(ax.get_xticklabels(), fontsize=8)	#rotation='vertical', 
			plt.setp(ax.get_yticklabels(), fontsize=8)
			plt.xlabel("Germline Divergence (%)", fontsize=8)
			plt.ylabel("SeqID to %s (%%)"%native, fontsize=8)
			
			#cbar = fig.colorbar(cax)


			cax 	= ax.imshow(my_array, vmax=4, vmin=0, origin='lower')
			plt.tight_layout()
			
			"""
			#plt.contourf(my_array)
			# here add native identities
			if native in nat_set and is_light == 0:
				#print "loading native identity/divergence..."
				nats, xs, ys = load_divid(native, is_light)
				ax.plot(xs, ys, "ko", marker='.', alpha=1, markersize=6)
				
				
				# add text to the plot
				xyts = sorted(zip(ys, xs, nats))
				for index, (y, x, t) in enumerate(xyts):
					shift = 3 * len(t)
					if t == native:
						ax.text(x + 1, y - 3.5, t, fontsize=8, fontweight="bold")
					else:
						if index % 2 == 0:
							ax.text(x + 0.6, y - 1.5, t, fontsize=8, fontweight="bold")
						else:
							ax.text(x - shift - 1.5, y - 1.5, t, fontsize=8, fontweight="bold")
				
			"""	
			
			
						
			#sys.exit(0)
			print "saving..."
			plt.savefig("%s/%s-%s_%s.png" %(prj_tree.figure, f_head, native, title), dpi=600)
			#plt.show()
			#print "done!"
			
			#break
			
			


if __name__ == '__main__':
	# get parameters from input
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, l="is_light", r="reads", t="title")
	is_light, reads, title = getParasWithDefaults(dict_args, dict(),"is_light","reads","title") #use an empty dict for defaults rather than rewrite general method to avoid errors

	r_list = list()
	if reads is not None:
		with open(reads) as x: r_list = map(str.strip, x.readlines())

	prj_tree 	= ProjectFolders(os.getcwd())
	prj_name 	 = fullpath2last_folder(prj_tree.home)
	#fig_title = DICT_DATASET_TITLE[prj_name]
	

	nat_set = set(["VRC01", "VRC03", "VRC06", "VRC08", "VRC07"])
	if is_light == 1:
		nat_set = set(["VRC01", "VRC03", "VRC06", "VRC07"])
	elif is_light >= 10:
		nat_set = set(["CH505-1", "CH505-2", "CH505-4"])
	
	main()

