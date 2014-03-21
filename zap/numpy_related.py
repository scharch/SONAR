#!/usr/bin/env python
# encoding: utf-8
from numpy import mean, array, zeros, ones, nan, std, isnan
from math import log

def generate_array(dict_div_id):
	my_array, total = zeros( (10 ** 2, 10 ** 2) ), 0

	for key, (divergence, identity) in dict_div_id.items():
		
		divergence 	= int(round(divergence, 0))
		identity 	= int(round(identity, 0))
		
		div_ind = int(divergence)
		id_ind 	= int(identity)
		
		try:
			my_array[id_ind, div_ind] += 1
			total += 1
		except:
			pass;			# no identity information
				
	return my_array	
	

def load_read2germ_assignment(align_file , germ):
	print "loading reads assigned to germline: %s from file: %s..." %(germ, align_file)
	result, total 	= dict(), 0
	reader			= csv.reader(open(align_file, "rU"), delimiter = sep)
	for row in reader:
		qid, sid = row[0].strip(), row[1].strip()
		if sid == germ:
			result[qid] = [nan] * 2
			total += 1

	print "%d reads loaded..." %total
	return result

def manipulate_array(my_array):
	
	x, y = my_array.shape
	for i in range(x):
		for j in range(y):
			if my_array[i, j] == 0.0:
				my_array[i, j] = nan
			elif my_array[i, j] == 1.0:
				my_array[i, j] = 0.01
			else:
				my_array[i, j] = log(my_array[i, j], 10)
	
	return my_array


