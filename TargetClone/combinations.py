#Handler for classes that involve making combinations

from alleles import Alleles
from gaussianMM import GaussianMM
import numpy as np
import itertools
import sys
import settings

#Define the combinations in $C_c$ that can be made. For every 2 measurement positions in 2 subclones, we have 4 positions in a matrix where we can vary the copy number between kmin and kmax.
#The copy number of the parent at measurement position i and j are always given (known from the previous iteration), but for the first measurement position the first two positions need to be inferred at once. We here make
#all possible combinations between kmin and kmax for these two positions in the subclone. 
class CCombinations:
	combinations = None
	
	def __init__(self, kmin, kmax):
		self.combinations = self.defineCombinations(kmin, kmax)
		
	def defineCombinations(self, kmin, kmax):
		#Make the combinations with C objects. The scores can then still all be pre-computed.
		#Version with only 2^6 combinations
		return list(itertools.product(range(kmin,kmax+1), range(kmin, kmax+1)))
	
	# #Given a laf we can compute the best combination. 
	# def getBestCombination(self, laf):
	# 	return None

#Handler for making combinaions between alleles. A combination between alleles is a combination between the alleles of the normal cell and the alleles of the tumor subclone. 	
class AlleleCombination:

	alleles = None #holds the normal alleles and the tumor alleles
	eventDistance = 0
	laf = None
	
	def __init__(self, alleles):
		self.alleles = alleles
	
	#Compute the event distance for this combination given the alleles of the parent.
	#We compute the event distance only for the tumor alleles to the given parent allele
	def computeEventDistance(self, parentAlleles):
		for allele in range(1, len(self.alleles)):
			self.eventDistance += self.alleles[allele].computeEditDistance(parentAlleles)
		
		return 1/float(self.eventDistance + 1) 
			
	#Compute the LAF of this combination of alleles given the frequency of the healthy and tumor subclone
	def computeLAF(self, mu):
		#Compute the LAF that this combination of alleles would generate given a mu in a mixture
		totalACounts = 0
		totalBCounts = 0
	
		for allele in range(0, len(self.alleles)): #THe first position corresponds to the normal subclone, the second position to the tumor subclone
			totalACounts += (self.alleles[allele].ACount * mu.mu[allele])
			totalBCounts += (self.alleles[allele].BCount * mu.mu[allele])
			
			
		totalCount = float(totalACounts + totalBCounts)	
		if totalCount == 0:
			self.laf = 0
		else:
			self.laf = min(totalACounts, totalBCounts) / totalCount

		return self.laf
		
#Handler for functions on combinations of C and Mu. This class is mainly used to compute the probability of observing a certain combination of C and mu given LAF measurements. 		
class CMuCombination:

	c = None
	mu = None
	mixtureModels = None #There are multiple mixture models, one for each possible parent
	probabilities = None #the probabilities for every possible LAF value (including noise) for this C and mu combination
	mixtureModelValleys = None #For each mixture model (for every parent), we pre-compute the valleys to save time. From these we can determine which alleles will be assigned to a LAF. 
	
	def __init__(self, c, mu, eventDistances):
		self.c = c
		self.mu = mu
		self.probabilities = []
		self.pLafGivenMu = []
		self.computeCombinationLaf()
		self.computeLafMixtureModels(eventDistances)
		self.setValleys()
		
	#Compute the LAF that we get from combining this C + mu	
	def computeCombinationLaf(self):
		for i in range(0, len(self.c.alleleCombinations)):
			self.c.alleleCombinations[i].computeLAF(self.mu) #each C can have multiple possible allele explanations, each of these has its own laf
			
	#Pre-compute the mixture model of probabilities that we would get for this C and mu combination. If this information is readily available, we only need to query it with a laf to get a score		
	def computeLafMixtureModels(self, eventDistances):
		self.mixtureModels = {}
		combinations = self.c.alleleCombinations
		
		#1. Build the mixture density from lambda, mu & sigma of the Alleles object
		sigma = settings.general['mmSigma'] #this will need to be estimated somewhere
		
		#2. The probabilities are first set at each LAF with the normalized event distance as a probability value
		for parentAllele in eventDistances.alleleList:
			
			#obtain the event distances to the combinations of this C
			totalDistances = []
			for combination in combinations:
				combinationAlleles = combination.alleles[1].alleleString
				if self.mu.mu[0] == 1: #if there are only normal cells, we compute the distance based on the normal cells. 
					combinationAlleles = 'AB'
				#obtain the normalized event distance for this set of alleles (normal combined with tumor)
				dist = eventDistances.getEventDistance(parentAllele, combinationAlleles)
			
				totalDistances.append(dist)
			#compute the ratios of the mixture model
			ratios = np.divide(totalDistances, np.sum(totalDistances))
		
			#define the mixture model
			self.mixtureModels[parentAllele] = GaussianMM(self.c.getCombinationLafs(), sigma, ratios)

	
	def addProbabilities(self, probability):
		self.probabilities.append(probability)
	
	#obtain the mixture model for a given parent. A different set of parental alleles will result in a different event distance and thus a different mixture model. Thus, we need to select the right model given the parent alleles	
	def getMixtureModelByParent(self, parentAlleles):
		return(self.mixtureModels[parentAlleles])


	#we need to go through each possible mixture model depending on the parent, and then pre-compute what the valleys will be. 
	def setValleys(self):
		self.mixtureModelValleys = dict()
		for parentMixtureModelKey in self.mixtureModels:
			
			#Obtain the mixture model in question
			mixtureModel = self.mixtureModels[parentMixtureModelKey]
			valleys = mixtureModel.getValleys() #these are the tresholds, they are assumed to be in the valleys. Actual valley detection is tricky with a lot of sequencing noise (individual distributions overlap)
		
			#Get the number of the valley that the laf is located in, use this to get the allele combination back
			startLaf = [0]
			endLaf = [0.5]
			
			moddedValleys = startLaf + valleys + endLaf
			
			self.mixtureModelValleys[parentMixtureModelKey] = moddedValleys


	###! This function is the one that slows down the code the most. Can we optimize it? Which parts are the slowest? 

	#obtain the alleles most likely corresponding to a given LAF, for this we use the allele threshold.
	#The threshold (or multiple) is defined at the mean between every pair of two possible LAF that we can measure with this C and Mu combination.
	#Depending on if the measured LAF is on the right or left of the treshold, the most likely set of alleles is reported. 
	def getAllelesByLaf(self, parentAlleles, laf): #is this the right place for this step?

		if self.mu.mu[0] == 1: #if there is no tumor component the alleles can only be AB by assumption
			return Alleles(1,1)
			
		#Old code to determine the valleys, this is now pre-computed. 
		#we need to know the parent alleles to get the right mixture model to extract the right alleles given a laf
		# mixtureModel = self.getMixtureModelByParent(parentAlleles.getAllelesAsString())
		# 
		# #if we know where the local minimum starts and ends, we know which regions are associated with which alleles
		# #In the event of no minima, the array of valleys will be empty.
		# 
		# #Getting the valleys is a very slow function, slowdown with ~ 2 seconds
		# valleys = mixtureModel.getValleys() #these are the tresholds, they are assumed to be in the valleys. Actual valley detection is tricky with a lot of sequencing noise (individual distributions overlap)
		# 
		# #Get the number of the valley that the laf is located in, use this to get the allele combination back
		# startLaf = [0]
		# endLaf = [0.5]
		# 
		# moddedValleys = startLaf + valleys + endLaf
		
		#here we can now just lookup the valleys without having to compute them
		moddedValleys = self.mixtureModelValleys[parentAlleles.getAllelesAsString()]
		
		#between which possible LAF is the measurement value located? THe alleles remain ordered, so we can select the right allele combination. 
		alleleInd = 0
		for border in range(0,len(moddedValleys)-1):
		  	if laf >= moddedValleys[border] and laf <= moddedValleys[border+1]:
		  		alleleInd = border
		
		return self.c.alleleCombinations[alleleInd].alleles[1]
		
		
		