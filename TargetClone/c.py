import numpy as np
from combinations import AlleleCombination
from alleles import Alleles

#Handler for operations on the C matrix. In this context, C always has one row, so every C class represents one measurement value and one normal sample and one tumor subclone (2x1, cxr). 
class C:

	c = None
	alleles = None #the total set of alleles
	distance = None
	alleleCombinations = None #every chromosome copy at a measurement position is assumed to correspond to one allele, so each measurement position has multiple possible combination of alleles given the amount of copies in C. 
	
	#Initialize the class for run speed
	def __init__(self, c = None, alleleCombinations = None):
		if c is None:
			self.c = []
		else:
			self.c = c
		if alleleCombinations is None:
			self.alleleCombinations = []
			self.setAlleleCombinations()
		else:
			self.alleleCombinations = []
	
	#set C given the initial value of the tumor subclone. The normal component in a sample always has copy number 2. 
	def setC(self, cInit):
		self.c[0] = 2
		self.c[1] = cInit
	
	#Given the copy number of the tumor subclone, we can obtain all the possible combinations of alleles that we could have observed
	def setAlleleCombinations(self):
		k = self.c[1]
		for copyN in range(0,k+1):
			ACount = copyN
			BCount = k - copyN
			tumorAlleles = Alleles(ACount, BCount)
			normalAlleles = Alleles(1,1)
			self.alleleCombinations.append(AlleleCombination([normalAlleles, tumorAlleles]))
	
	#For every combination of alleles that is possible, the event distance to the normal sample is computed	
	def getCombinationEventDistances(self, parentAlleles):
		eventDistances = []
		for i in range(0, len(self.alleleCombinations)):
			eventDistances.append(self.alleleCombinations[i].computeEventDistance(parentAlleles))
		
		return eventDistances
	
	#The LAF is computed for every possible combination of alleles that we can make with a copy number
	def getCombinationLafs(self):
		lafs = []
		for i in range(0, len(self.alleleCombinations)):
			lafs.append(self.alleleCombinations[i].laf)

		return lafs
		