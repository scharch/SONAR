
"""

defines a GSSP class with methods for reading from/writing to text files
and calculating rarity, Jensen-Shannon divergence between profiles, and 
position-wise Shannon entropy

"""

from .. import *
import csv
from collections import defaultdict
from math import log
from numpy import average, std
import pandas

class GSSP:
	
	# spectrum (dict) -> vgenes (list) -> sample[s] (list)	-> positions (dict) -> freq (num), profile (list)

	def __init__(self, inFile, name="unnamed"):

		self.name    = name
		self.vgenes  = defaultdict( list )
		self.average = defaultdict( list )
		self.rarity  = []
		self.entropy = []

		#read in GSSPs from a text file
		with open(inFile, "rU") as handle:
			reader = csv.reader(handle, delimiter = "\t")
			header = next(reader)
			thisV  = None
			sample = []
			for row in reader:
				if int(row[2]) == 1:
					if thisV is not None:
						self.vgenes[thisV].append( sample )
					sample = []
					thisV = row[0]
				sample.append( dict( germline=row[3].split(","), freq=eval(row[4]), profile=list(map(float, row[5:])) ) )
			if thisV is not None: #relevant for empty profile (no genes with enough seqs)
				self.vgenes[thisV].append( sample )


	#this function compares the GSSPs from all V genes between two sets of profiles
	def compare(self, otherSpectrum, subset=""):
		results = defaultdict( dict )
		for v1 in self.vgenes:
			if subset not in v1:
				continue
			for v2 in otherSpectrum.vgenes:
				if subset not in v2:
					continue
				comps = []
				for i in self.vgenes[v1]:
					for j in otherSpectrum.vgenes[v2]:
						comps.append( spectrumJSD(i,j,v1,v2) )
				results["%s,%s"%(self.name,v1)]["%s,%s"%(otherSpectrum.name,v2)] = average(comps)

		return pandas.DataFrame(results)
		
	
	#this function compares GSSPs between V genes in the current set of profiles
	def betweenV(self, subset=""):
		results = defaultdict( dict )
		k = sorted( self.vgenes.keys() )
		for n, v1 in enumerate(k):
			if subset not in v1:
				continue
			results["%s,%s"%(self.name,v1)]["%s,%s"%(self.name,v1)] = 0 #by definition
			for v2 in k[n+1:]:
				if subset not in v2:
					continue
				comps = []
				for i in self.vgenes[v1]:
					for j in self.vgenes[v2]:
						comps.append( spectrumJSD(i,j,v1,v2) )
				results["%s,%s"%(self.name,v1)]["%s,%s"%(self.name,v2)] = average(comps)
				results["%s,%s"%(self.name,v2)]["%s,%s"%(self.name,v1)] = average(comps)

		return pandas.DataFrame(results)


	def computeRarity(self, fixFreq=1):

		#the fixFreq argument is because I screwed up in the calculations of frequency
		#for spectrums generated before 2016-06-02
		#it's a constant factor, at least so only matters for rarity
		#can remove this once I re-generate data...

		aaList = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		for v in sorted(self.vgenes.keys()):
			r = defaultdict( lambda: defaultdict( list ) )
			for i in self.vgenes[v]:
				for p, pos in enumerate(i):
					if pos['freq'] is None:
						continue
					for aa, mut in enumerate(pos['profile']):
						if aaList[aa] in pos['germline']:
							continue #skip germline residues
						r[p][aaList[aa]].append( 1 - (pos['freq']*mut/fixFreq) )
			#go back through it so we can get position-wise averages over samples
			for p in r:
				for aa in r[p]:
					self.rarity.append( [ self.name, v, p+1, ",".join(i[p]['germline']), aa, "%.3f"%average(r[p][aa]), "%.3f"%std(r[p][aa]) ] )


	#averages multiple GSSPs sampled from a single gene
	# currently 
	def averageProfile(self):
		for v in self.vgenes:
			r = [ [0] * 20 for dummy in range(len(self.vgenes[v][0])) ]
			f = [ []  for dummy in range(len(self.vgenes[v][0])) ]
			for i in self.vgenes[v]:
				for p, pos in enumerate(i):
					if pos['freq'] is None:
						continue
					r[p] = [ r[p][a] + aa for a,aa in enumerate(pos['profile']) ]
					f[p].append(pos['freq'])
			#now divide out and store
			for x in range(len(f)):
				if len(f[x]) == 0:
					#masked position
					self.average[v].append( dict( freq=None, profile=[ 0.00 for y in range(20) ] ) )
				else:
					self.average[v].append( dict( freq=sum(f[x])/len(f[x]), profile=[ r[x][y]/len(f[x]) for y in range(20) ] ) )


	def profileEntropy(self, use_all=True):
		if use_all:
			if len(self.average) == 0:
				self.averageProfile()
			for v in sorted(self.average.keys()):
				e = []
				w = []
				for pos in self.average[v]:
					if pos['freq'] is None:
						continue
					e.append( shannon(pos['profile']) )
					w.append( pos['freq'] )
				self.entropy.append( [ self.name, v, "%.3f"%average(e, weights=w) ] )
		else:
			for v in sorted(self.vgenes.keys()):
				e = []
				w = []
				for pos in self.vgenes[v][0]:
					if pos['freq'] is None:
						continue
					e.append( shannon(pos['profile']) )
					w.append( pos['freq'] )
				self.entropy.append( [ self.name, v, "%.3f"%average(e, weights=w) ] )




def shannon( profile ):
	return -1 * sum([ float(x) * log(float(x),2) if x>0 else 0 for x in profile ])




def positionJSD( profile1, profile2, renorm=True ):

	if isinstance(profile2, str):
		profile2 = letter2profile(profile2)

	if renorm:
		profile1 = [ float(a) / sum(profile1) if sum(profile1)>0 else 1.0/len(profile1) for a in profile1 ]
		profile2 = [ float(a) / sum(profile2) if sum(profile2)>0 else 1.0/len(profile2) for a in profile2 ]

	averaged_profile = [ (float(a)+float(b))/2 for a,b in zip(profile1, profile2) ]

	return shannon(averaged_profile) - (shannon(profile1) + shannon(profile2))/2




def letter2profile( letter ):
	aa_order = dict( A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9, M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19 )
	profile = [0] * 20
	if letter != "X": profile[ aa_order[letter] ] = 1
	return profile




def spectrumJSD( spectrum1, spectrum2, v1=None, v2=None ):

	if not v1 == v2:
		#note - this dictionary is in position (1-based) counting...
		indels = { 'IGHV1-18'	:  dict( ins=[], dels=[] ),
			   'IGHV1-2'	:  dict( ins=[], dels=[] ),
			   'IGHV1-24'	:  dict( ins=[], dels=[] ),
			   'IGHV1-3'	:  dict( ins=[], dels=[] ),
			   'IGHV1-45'	:  dict( ins=[], dels=[] ),
			   'IGHV1-46'	:  dict( ins=[], dels=[] ),
			   'IGHV1-58'	:  dict( ins=[], dels=[] ),
			   'IGHV1-69'	:  dict( ins=[], dels=[] ),
			   'IGHV1-8'	:  dict( ins=[], dels=[] ),
			   'IGHV1-C'	:  dict( ins=[], dels=[] ),
			   'IGHV1-F'	:  dict( ins=[], dels=[] ),
			   'IGHV2-26'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV2-5'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV2-70'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV3-11'	:  dict( ins=[], dels=[] ),
			   'IGHV3-13'	:  dict( ins=[], dels=[54] ),
			   'IGHV3-15'	:  dict( ins=[54,55], dels=[] ),
			   'IGHV3-16'	:  dict( ins=[], dels=[] ),
			   'IGHV3-20'	:  dict( ins=[], dels=[] ),
			   'IGHV3-21'	:  dict( ins=[], dels=[] ),
			   'IGHV3-23'	:  dict( ins=[], dels=[] ),
			   'IGHV3-30'	:  dict( ins=[], dels=[] ),
			   'IGHV3-30-3' :  dict( ins=[], dels=[] ),
			   'IGHV3-33'	:  dict( ins=[], dels=[] ),
			   'IGHV3-35'	:  dict( ins=[], dels=[] ),
			   'IGHV3-38'	:  dict( ins=[], dels=[53,54] ),
			   'IGHV3-43'	:  dict( ins=[], dels=[] ),
			   'IGHV3-48'	:  dict( ins=[], dels=[] ),
			   'IGHV3-49'	:  dict( ins=[54,55], dels=[] ),
			   'IGHV3-53'	:  dict( ins=[], dels=[54] ),
			   'IGHV3-64'	:  dict( ins=[], dels=[] ),
			   'IGHV3-66'	:  dict( ins=[], dels=[54] ),
			   'IGHV3-7'	:  dict( ins=[], dels=[] ),
			   'IGHV3-72'	:  dict( ins=[59,60,61], dels=[54] ),
			   'IGHV3-73'	:  dict( ins=[59,60,61], dels=[54] ),
			   'IGHV3-74'	:  dict( ins=[], dels=[] ),
			   'IGHV3-9'	:  dict( ins=[], dels=[] ),
			   'IGHV3-D'	:  dict( ins=[], dels=[53,54] ),
			   'IGHV3-NL1'	:  dict( ins=[], dels=[] ),
			   'IGHV4-28'	:  dict( ins=[31], dels=[54] ),
			   'IGHV4-30-2' :  dict( ins=[31,32], dels=[54] ),
			   'IGHV4-30-4' :  dict( ins=[31,32], dels=[54] ),
			   'IGHV4-31'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV4-34'	:  dict( ins=[], dels=[54] ),
			   'IGHV4-39'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV4-4'	:  dict( ins=[31], dels=[54] ),
			   'IGHV4-59'	:  dict( ins=[], dels=[54] ),
			   'IGHV4-61'	:  dict( ins=[31,32], dels=[54] ),
			   'IGHV4-B'	:  dict( ins=[31], dels=[54] ),
			   'IGHV5-51'	:  dict( ins=[], dels=[] ),
			   'IGHV5-A'	:  dict( ins=[], dels=[] ),
			   'IGHV6-1'	:  dict( ins=[31,32,61,62], dels=[54] ),
			   'IGHV7-4-1'	:  dict( ins=[], dels=[] ),
			   'IGHV7-81'	:  dict( ins=[], dels=[] ),
			   'IGKV1-12'	:  dict( ins=[], dels=[] ),
			   'IGKV1-16'	:  dict( ins=[], dels=[] ),
			   'IGKV1-17'	:  dict( ins=[], dels=[] ),
			   'IGKV1-27'	:  dict( ins=[], dels=[] ),
			   'IGKV1-33'	:  dict( ins=[], dels=[] ),
			   'IGKV1-37'	:  dict( ins=[], dels=[] ),
			   'IGKV1-39'	:  dict( ins=[], dels=[] ),
			   'IGKV1-5'	:  dict( ins=[], dels=[] ),
			   'IGKV1-6'	:  dict( ins=[], dels=[] ),
			   'IGKV1-8'	:  dict( ins=[], dels=[] ),
			   'IGKV1-9'	:  dict( ins=[], dels=[] ),
			   'IGKV2-24'	:  dict( ins=[30,31,32,33,34], dels=[] ),
			   'IGKV2-28'	:  dict( ins=[30,31,32,33,34], dels=[] ),
			   'IGKV2-30'	:  dict( ins=[30,31,32,33,34], dels=[] ),
			   'IGKV2-40'	:  dict( ins=[30,31,32,33,34,35], dels=[] ),
			   'IGKV3-11'	:  dict( ins=[], dels=[] ),
			   'IGKV3-15'	:  dict( ins=[], dels=[] ),
			   'IGKV3-20'	:  dict( ins=[30], dels=[] ),
			   'IGKV3-7'	:  dict( ins=[30], dels=[] ),
			   'IGKV4-1'	:  dict( ins=[30,31,32,33,34,35], dels=[] ),
			   'IGKV5-2'	:  dict( ins=[], dels=[] ),
			   'IGKV6-21'	:  dict( ins=[], dels=[] ),
			   'IGKV6-41'	:  dict( ins=[], dels=[] ),
			   'IGLV1-36'	:  dict( ins=[], dels=[] ),
			   'IGLV1-40'	:  dict( ins=[33], dels=[] ),
			   'IGLV1-41'	:  dict( ins=[], dels=[] ),
			   'IGLV1-44'	:  dict( ins=[], dels=[] ),
			   'IGLV1-47'	:  dict( ins=[], dels=[] ),
			   'IGLV1-50'	:  dict( ins=[33], dels=[] ),
			   'IGLV1-51'	:  dict( ins=[], dels=[] ),
			   'IGLV10-54'	:  dict( ins=[], dels=[] ),
			   'IGLV11-55'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV2-11'	:  dict( ins=[33], dels=[] ),
			   'IGLV2-14'	:  dict( ins=[33], dels=[] ),
			   'IGLV2-18'	:  dict( ins=[33], dels=[] ),
			   'IGLV2-23'	:  dict( ins=[33], dels=[] ),
			   'IGLV2-33'	:  dict( ins=[33], dels=[] ),
			   'IGLV2-8'	:  dict( ins=[33], dels=[] ),
			   'IGLV3-1'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-10'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-12'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-16'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-19'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-21'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-22'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-25'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-27'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-32'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV3-9'	:  dict( ins=[], dels=[27,28] ),
			   'IGLV4-3'	:  dict( ins=[31,54,55,56,57], dels=[31,32] ),
			   'IGLV4-60'	:  dict( ins=[31,54,55,56,57], dels=[31,32] ),
			   'IGLV4-69'	:  dict( ins=[31,54,55,56,57], dels=[31,32] ),
			   'IGLV5-37'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV5-39'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV5-45'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV5-48'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV5-52'	:  dict( ins=[33,56,57,58,59,73,74], dels=[] ),
			   'IGLV6-57'	:  dict( ins=[68,69], dels=[] ),
			   'IGLV7-43'	:  dict( ins=[33], dels=[] ),
			   'IGLV7-46'	:  dict( ins=[33], dels=[] ),
			   'IGLV8-61'	:  dict( ins=[33], dels=[] ),
			   'IGLV9-49'	:  dict( ins=[31], dels=[31,32] )
			       }
		
		try:
			#kludgy but: use set notation to get differences (so if comparing 4-4 to 4-28, then no need to change anything)
			#	     then sort max to min so we don't shift indices when splicing
			#	     For the same reason, kills ins before accounting for dels

			ins1 = indels[v1]['ins']
			ins2 = indels[v2]['ins']
			for i in sorted( set(ins1).difference(ins2), reverse=True ):
				spectrum1 = spectrum1[0:i-1] + spectrum1[i:]
			for i in sorted( set(ins2).difference(ins1), reverse=True ):
				spectrum2 = spectrum2[0:i-1] + spectrum2[i:]

			dels1 = indels[v1]['dels']
			dels2 = indels[v2]['dels']
			for d in sorted( set(dels1).difference(dels2), reverse=True ):
				spectrum2 = spectrum2[0:d-1] + spectrum2[d:]
			for d in sorted( set(dels2).difference(dels1), reverse=True ):
				spectrum1 = spectrum1[0:d-1] + spectrum1[d:]

		except KeyError:
			print( "Warning: Unrecognized V gene attempting to compare %s and %s - YOUR RESULTS WILL BE WRONG!!" % (v1, v2) )


	currentJSD = []
	weights	   = []
	for p1, p2 in zip(spectrum1, spectrum2):
		if p1['freq'] is not None and p2['freq'] is not None:
			currentJSD.append( positionJSD( p1['profile'], p2['profile'] ) )
			weights.append( average( [ p1['freq'], p2['freq'] ] ) )

	return average( currentJSD, weights=weights )
