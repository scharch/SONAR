
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
		self.rarity  = defaultdict( dict )
		self.entropy = dict()

		#read in GSSPs from a text file
		with open(inFile, "r") as handle:
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
	def compare(self, otherSpectrum, alignment="<SONAR>/sample_data/GSSPs/comparison_dictionary.csv", subset=""):
		alignment = re.sub( "<SONAR>", SCRIPT_FOLDER, alignment)
		alignDict = dict()
		with open(alignment, 'r') as handle:
			alignCSV = csv.reader(handle)
			for row in alignCSV:
				alignDict[ row[0] ] = dict( ins=[int(x) for x in row[1].split(",") if not row[1]==""], dels=[int(x) for x in row[2].split(",") if not row[2]==""] )
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
						comps.append( spectrumJSD(i,j,alignDict,v1,v2) )
				results["%s,%s"%(self.name,v1)]["%s,%s"%(otherSpectrum.name,v2)] = average(comps)

		return pandas.DataFrame(results)
		
	
	#this function compares GSSPs between V genes in the current set of profiles
	def betweenV(self, alignment="<SONAR>/sample_data/GSSPs/comparison_dictionary.csv", subset=""):
		alignment = re.sub( "<SONAR>", SCRIPT_FOLDER, alignment)
		alignDict = dict()
		with open(alignment, 'r') as handle:
			alignCSV = csv.reader(handle)
			for row in alignCSV:
				alignDict[ row[0] ] = dict( ins=[int(x) for x in row[1].split(",") if not row[1]==""], dels=[int(x) for x in row[2].split(",") if not row[2]==""] )
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
						comps.append( spectrumJSD(i,j,alignDict,v1,v2) )
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
				self.rarity[ v ][ p ] = dict( germline=",".join(i[p]['germline']), mutants=dict() )
				for aa in r[p]:
					self.rarity[ v ][ p ][ 'mutants' ][ aa ] = dict( average=average(r[p][aa]), stddev=std(r[p][aa]) )


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
				self.entropy[ v ] = average(e, weights=w) 
		else:
			for v in sorted(self.vgenes.keys()):
				e = []
				w = []
				for pos in self.vgenes[v][0]:
					if pos['freq'] is None:
						continue
					e.append( shannon(pos['profile']) )
					w.append( pos['freq'] )
				self.entropy[ v ] = average(e, weights=w) 



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




def spectrumJSD( spectrum1, spectrum2, indels, v1=None, v2=None ):

	if not v1 == v2:
		#note - this dictionary is in position (1-based) counting...
		
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
