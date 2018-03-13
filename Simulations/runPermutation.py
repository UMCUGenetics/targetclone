import sys
sys.path.insert(0, '../TargetClone/')

from multiprocessing import Pool
from copy import deepcopy
import random
import numpy as np
import pickle
import simulationSettings
import sys
from laf import LAF
from mu import Mu

from simulationErrors import SimulationErrorHandler

class Permute:
	
	positions = None
	chromosomes = None
	
	simulationData = None
	samples = None
	simulator = None
	targetClone = None
	savedMu = None #store the mu separately, apparently these are somehow pdated by TC even though we use deepcopy
	
	cScores = None
	aScores = None
	muScores = None
	treeScores = None
	
	pC = None
	pA = None
	pMu = None
	pTrees = None
	pSamples = None
	
	uniqueIdentifier = None
	
	#We need: the AF matrix and SNV matrix to permute
	#The real C, 
	
	#2017-05-16T11:08:22.753351
	def __init__(self):
		#print "reached"
		uniqueIdentifier = sys.argv[1]
		print uniqueIdentifier
		self.uniqueIdentifier = uniqueIdentifier
		
		
		#Obtain: data locations in a settings file, file names by unique identifiers
		#Parse the data from these locations
		
		self.positions = []
		self.chromosomes = []
		
		self.cScores = []
		self.aScores = []
		self.muScores = []
		self.treeScores = []
		
		self.pC = []
		self.pA = []
		self.pSamples = []
		self.pTrees = []
		self.pMu = []
		
		#instantiate TC object
		targetCloneInstance = simulationSettings.files['targetCloneInstance']
		pkl_file = open(targetCloneInstance, 'rb')
		self.targetClone = pickle.load(pkl_file)
		pkl_file.close()
		
		#1. initialize

		self.parseSimulationData(uniqueIdentifier)
		
		#2. Run a permutation
		[pCMatrix, pAMatrix, pSamples, pTreeList, permutedC, permutedA] = self.permute()
		#3. Compute the errors
		[pCError, pAError, pMuError, pTreeScore, ambiguityScore, totalScore] = self.computePermutationErrors(pCMatrix, pAMatrix, pSamples, pTreeList, permutedC, permutedA)
		#4. Write the errors to a file
		self.writeErrorsToFiles(uniqueIdentifier, pCError, pAError, pMuError, pTreeScore, ambiguityScore, totalScore)
		
		
		
	def readDataFromFile(self, file):
		text_file = open(file, "r")
		lines = text_file.read()
		floatLines = []
		
		for line in lines.split("\n"):
			if line != "":
				floatLines.append(float(line))
		
		text_file.close()
		
		return floatLines
	
	#It is not efficient to re-do this every time. Maybe do it once and then use pickle to store the matrices? 
	def parseSimulationData(self, uniqueIdentifier):
		#Load the file by unique identifier, and store the object
		#Get the actual location from the settings file!
		pkl_file = open(simulationSettings.files['outputDir'] + uniqueIdentifier + '/simulationData.pkl', 'rb')
		self.simulationData = pickle.load(pkl_file)
		pkl_file.close()
		
		print self.simulationData.realTree.edgeList
		
		
		#We will have to obtain the real mu from the file.
		muFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/RealMu_0.txt'
		self.savedMu = self.readDataFromFile(muFile)
	
				
	def generateLAFObject(self, measurements):
		return LAF(measurements, self.chromosomes, self.positions, self.positions)
	
	def unison_shuffled_copies(self, a, b):
		assert len(a) == len(b)
		p = np.random.permutation(len(a))
		return a[p], b[p]

	def permute(self):
		
		#1. Shuffle A and AF in unison
		
		permutedA = deepcopy(self.simulationData.aMatrix)
		#to permute A, we need to use the indices since it is not working on objects!!
		indices = np.arange(0, permutedA.size)
		reshapedIndices = indices.reshape(permutedA.shape)
		
		# [permutedA1, permutedAF1] = self.unison_shuffled_copies(permutedA, self.simulationData.afMatrix)
		# [permutedA2, permutedAF2] = self.unison_shuffled_copies(np.transpose(permutedA1), np.transpose(permutedAF1))
		
		#original code
		#[permutedAF1, reshapedIndices1] = self.unison_shuffled_copies(self.simulationData.afMatrix, reshapedIndices)
		#[permutedAF2, reshapedIndices2] = self.unison_shuffled_copies(np.transpose(permutedAF1), np.transpose(reshapedIndices1))
		
		[permutedAF2, reshapedIndices2] = self.unison_shuffled_copies(self.simulationData.afMatrix, reshapedIndices)
		
		#permutedA2 = np.transpose(permutedA2)
		#First do an ind2sub to get back the allele information
		reshapedIndices2 = reshapedIndices2
		row,col = np.unravel_index(reshapedIndices2, reshapedIndices2.shape)
		
		permutedA2 = permutedA[row,col]
		
	
		#we need to permute the other matrices in the same way, so that the order between them corresponds. 
		
		permutedC = np.zeros(permutedA2.shape, dtype=int) #we can easily derive C from A without making a copy
		aString = np.empty(permutedA2.shape, dtype=object)
		#2. Derive C and AF from the shuffled A
		for row in range(0, permutedA2.shape[0]):
			for col in range(0, permutedA2.shape[1]):
				currentAllele = permutedA2[row][col].getAllelesAsString()
				#Count the length of the string
				currentAlleleC = len(currentAllele)
				permutedC[row][col] = currentAlleleC
				aString[row][col] = currentAllele
				
				#obtain the mu of the sample
				#use the mu to compute the AF measurement
				
		#write C and A to a file to test
		np.savetxt(simulationSettings.files['outputDir'] + self.uniqueIdentifier + '/PermutedC.txt', permutedC.astype(int), fmt='%i', delimiter='\t')
		
		#Save A
		np.savetxt(simulationSettings.files['outputDir'] + self.uniqueIdentifier + '/PermutedA.txt', aString, fmt='%s', delimiter='\t')

		#Save the permuted AF
		np.savetxt(simulationSettings.files['outputDir'] + self.uniqueIdentifier + '/PermutedAF.txt', permutedAF2, fmt='%s', delimiter='\t')
		
		#No SNVs, then the error for trees is very high in the normal simulation (9), same in the permutation. 
		#Removing edges does increase error
		#SNVs perfectly match what they should
		#Random SNVs decrease performance by a lot. 
		
		#I think it may be wise to do a test where SNVs do not influence the results. 
		
		permutedSamples = []
		
		#Convert the AF back to a LAF object
		#Set the somatic variants and other properties of the permuted samples. 
		for sampleInd in range(0, len(self.simulationData.samples)):
			permutedSample = deepcopy(self.simulationData.samples[sampleInd])
			permutedSample.afMeasurements = permutedAF2[:,sampleInd]
			
			newLAF = self.generateLAFObject(permutedAF2[:,sampleInd])
			permutedSample.measurements = newLAF
			permutedSamples.append(permutedSample)
		
		snvNum = len(permutedSamples[0].somaticVariants)
		snvMatrix = np.zeros([snvNum, len(permutedSamples)], dtype=int)
		
		from random import randint
		
		for sample in range(0, len(permutedSamples)):
			
			
			somVar = permutedSamples[sample].somaticVariants
			
			snvMatrix[:,sample] = somVar
			

		np.savetxt(simulationSettings.files['outputDir'] + self.uniqueIdentifier + '/PermutedSNV.txt', snvMatrix, fmt='%s', delimiter='\t')
		
		[pCMatrix, pAMatrix, pSamples, pTreeList, iterationMu, iterationMessages] = self.targetClone.run(permutedSamples)
		print "edges: "
		print pTreeList
		
		for tree in pTreeList:
			print pTreeList[tree].edgeList
		#Also return the original permuted matrices so that we can compute the error compared to these matrices!
		return [pCMatrix, pAMatrix, pSamples, pTreeList, permutedC, permutedA2]
	
	def generateLAFObject(self, afMeasurements):
		measurements = []
		for measurement in afMeasurements:
			if measurement > 0.5:
				measurements.append(1-measurement)
			else:
				measurements.append(measurement)
		return LAF(measurements, self.simulationData.chromosomes, self.simulationData.positions, self.simulationData.positions)
		
	
	def computePermutationErrors(self, pCMatrix, pAMatrix, pSamples, pTreeList, permutedC, permutedA):
		
		#Import a class specific for computing permutation errors
		#Use it to compute the errors below
		simulationErrorHandler = SimulationErrorHandler()
		
		#Obtain the muData from the samples (also provide the real mu!)
		realSamples = self.simulationData.samples
		pMuTumor = []
		realMuTumor = []
		
		for sample in range(0, len(pSamples)):
			pSampleMuT = pSamples[sample].bestCMu[0].mu.mu[1]
			#sampleMuT = realSamples[sample].bestCMu[0].mu.mu[1]
			
			#print "current saved mu: ", self.savedMu[sample]
			#print "normal mu: ", 1-(100*self.savedMu[sample])
			tumorMuConverted = int(100*self.savedMu[sample])
			realMuTumor.append(Mu(100-tumorMuConverted))
			#sampleMuT = realSamples[sample].bestCMu[0].mu
			pMuTumor.append(pSampleMuT)
			#savedMu.append(sampleMuT)
			
		
		[ambiguityScore, totalScore] = simulationErrorHandler.computeAmbiguityScore(permutedA, pAMatrix, self.savedMu, pMuTumor)
		
		pCMatrixFloat = pCMatrix.astype(float)
		error = simulationErrorHandler.computeCError(permutedC, pCMatrixFloat)
		pCError = error / float(permutedC.size)
		
		
		aData = simulationErrorHandler.computeAError(permutedA, pAMatrix)
		pAError = aData[0] / float(permutedA.size)
		
	
		muData = simulationErrorHandler.computeMuError(pSamples, realMuTumor)
		pMuError = muData[0]
		
		print self.simulationData.realTree.edgeList
		pTreeScore = simulationErrorHandler.computeTreeError(pTreeList, self.simulationData.realTree)
		
		print pCError
		print pAError
		print pMuError
		print pTreeScore
		print ambiguityScore
		print totalScore
		
		#return the scores
		return [pCError, pAError, pMuError, pTreeScore, ambiguityScore, totalScore]
	
	def writeErrorsToFiles(self, uniqueIdentifier, pCError, pAError, pMuError, pTreeScore, ambiguityScore, totalScore):
		#Get the path as specified in the settings
		#Write the results to these files (by unique identifier, append)
		pCFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pCError.txt'
		pAFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pAError.txt'
		pMuFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pMuError.txt'
		pTreeFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pTreeError.txt'
		pAmbiguityFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pAmbiguityError.txt'
		pAmbiguityCorrectedFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/pAmbiguityCorrectedError.txt'
		
		f = open(pCFile, 'a')
		f.write(str(pCError))
		f.write("\n")
		f.close()
		
		f = open(pAFile, 'a')
		f.write(str(pAError))
		f.write("\n")
		f.close()
		
		f = open(pMuFile, 'a')
		f.write(str(pMuError))
		f.write("\n")
		f.close()
		
		f = open(pTreeFile, 'a')
		f.write(str(pTreeScore))
		f.write("\n")
		f.close()
		
		f = open(pAmbiguityFile, 'a')
		f.write(str(ambiguityScore))
		f.write("\n")
		f.close()
		
		f = open(pAmbiguityCorrectedFile, 'a')
		f.write(str(totalScore))
		f.write("\n")
		f.close()
			
	
Permute()