import numpy as np
import time

#Handler for computing P(C,mu|LAF,T)
class Fitness:

	fitness = 0
	
	def computePCMuGivenLAF(self, cMuCombination, laf, parentAlleles):
		#1. Compute the P(C|T)
		#pC = self.computePC(cMuCombination.c)
		
		#2. Compute the P(LAF|C,mu,T)

		pLAFGivenCMu = self.computePLAFGivenCMu(cMuCombination,laf, parentAlleles)

		#3. Multiply the score
		#self.fitness = pC * pLAFGivenCMu #This is here turned off because we divide the distance by the distance for the whole combination of 4 positions in C_c, rather than for 1! 
		return pLAFGivenCMu
	
	def computePLAFGivenCMu(self, cMuCombination, laf, parentAlleles):
		#3. Obtain score from the pdf
		#Get the right mixture model for this parent
		return cMuCombination.mixtureModels[parentAlleles.alleleString].getPdf(laf)
	
	#Assume that this function is always provided a C object that is one row of C
	def computePC(self, C):
	
		# Get the distance for the current row (assume always the last row in C)
		dist = abs(C.c[0] - C.c[1])	+1	
		
		return (1/float(dist))
	