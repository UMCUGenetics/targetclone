"""
	Code to run a generic simulation case for TargetClone that is not biased towards TGCC. 

"""

##Imports

import sys
sys.path.insert(0, '../../TargetClone/')

import numpy as np
import csv
import re
import random
import uuid
#import matplotlib.pyplot as plt
#from pylab import axes
import scipy.stats as stats
from collections import OrderedDict
import math
import datetime
import os
import time
from multiprocessing import Pool
import pickle
from random import randint

from somaticVariants import SomaticVariant
from c import C
from alleles import Alleles
from copy import deepcopy
from sample import Sample
from mu import Mu
from laf import LAF
from combinations import CMuCombination
from run import TargetClone
from distance import EventDistances
from segmentations import Segmentation
from tree import Graph
from subclone import Subclone
from simulationDataParser import SimulationDataHandler
from simulationErrors import SimulationErrorHandler
from randomVectors import randvec
import settings
import simulationSettings

#Define intial parameters for this test, later replace with settings

#1. Do a subclonal expansion
class Simulator:
	
	snpNum = simulationSettings.general['numberOfSNPs']
	snvNum = simulationSettings.general['numberOfSNVs']
	simulationProbabilityFile = simulationSettings.files['simulationProbabilityFile']
	chromosomeArms = []
	allChromosomeArms = []
	simulationProbabilities = None
	snvs = [] #positions of all snv measurements
	kmin = simulationSettings.general['kmin']
	kmax = simulationSettings.general['kmax']
	muN = 0
	cellCycles = simulationSettings.general['cellCycles']
	usedVariantSlots = [] #Make sure that we do not use the same variant again in another subclone during evolution, store the used positions here. 
	uniqueID = str(randint(0,1000)) #will be overwritten by provided unique ID
	
	chromosomes = [] #All chromosomes of the SNP measurements (would have been better if this was an object too)
	positions = []
	
	
	def initializeGenericSimulation(self, muN, uniqueID): #allow for more parameters, then have switches for the different run types
		self.muN = muN
		self.uniqueID = uniqueID
		
		#1. Read the probabilities of losing/gaining chromosomes (or arms), these will be equal for the generic simulation but can be changed accordingly
		self.simulationProbabilities = SimulationProbabilities()
		
		self.simulationProbabilities.readFullProbabilityFile(self.simulationProbabilityFile)
		
		for armProbabilityList in self.simulationProbabilities.armProbabilities:
			self.chromosomeArms.append(armProbabilityList[0])
		
		#1.1 Distribute the SNP measurement positions across the chromosome arms (first equally)
		self.distributeSnpMeasurementsAcrossChromosomeArms()
		
		#1.2 Distirbute SNV positions randomly across the genome
		self.snvs = self.distributeSnvMeasurementsAcrossGenome()
		self.snvPositions = self.obtainSomaticVariantIndices()
		
		self.segmentation = Segmentation()
		self.segmentation.setSegmentationFromFile(simulationSettings.files['segmentationFile'])
		
		
		#2. Subclonal expansion
		
		self.subclonalExpansion()
		
	def distributeSnpMeasurementsAcrossChromosomeArms(self):
		
		#We do the distribution of SNP measurements when the number of desired SNPs is lower than the allowed boundary.
		#Otherwise, we will need to do a segmentation after that. Each chromosome arm will then have one SNP max. 
			
		#Distribute the SNPs across the genome based on the size of the chromosome. 
		#depending on the size of this chromosome. So Chromosome 1 has more than 22.
			
		#1. Get the chromosome sizes
		hg19CoordinateArray = self.getHg19Coordinates()
			
		
		
		#If the number of SNP measurements is too large, we recommend doing a segmentation. So that would mean that we have one measurement per chromosome arm,
		#assuming that the segmentation is perfect. Thus, we should distribute the measurements as such and update the snpNum value, as the rest of the code depends on that.
		#It should equal the number of chromosome arms.
		if simulationSettings.runType['segmentation'] == True:
			self.allChromosomeArms = []
			for armInd in range(0, len(self.chromosomeArms)):
				self.allChromosomeArms.append(self.chromosomeArms[armInd])
				
			self.snpNum = len(self.allChromosomeArms)
		else:
			
			
			sizes = []
			for arm in range(0, hg19CoordinateArray.shape[0]):
				
				size = hg19CoordinateArray[arm, 2] - hg19CoordinateArray[arm, 1]
				sizes.append(size)
			
			sizes = np.array(sizes)
			#2. Divide the size by the number of arms. Then multiply by the number of SNPs.
			#If this number is not an integer, round to the nearest integer
			 
			
			
			dividedSizes = sizes / float(np.sum(sizes))
			
			#Divide the SNPs across the chromosome arms based on their sizes (this should equal to the number of SNPs)
			snpDivision = np.round(dividedSizes * self.snpNum)
			
			self.allChromosomeArms = []
			emptyIndices = []
			for armInd in range(0, len(self.chromosomeArms)):
				arm = self.chromosomeArms[armInd]
				
				#Sometimes the addition of SNPs is more than we should add due to rounding
			
			#	print "previous length: ", len(self.allChromosomeArms)
				if len(self.allChromosomeArms) + int(snpDivision[armInd]) > self.snpNum:
					#print "new length: ", len(self.allChromosomeArms) + int(snpDivision[armInd])
					maxAddition = self.snpNum - len(self.allChromosomeArms)
					#print "max to add: ", maxAddition
					self.allChromosomeArms += [arm]*maxAddition
				else:
					self.allChromosomeArms += [arm]*int(snpDivision[armInd])
				
				#Sometimes not all chromosomes get SNPs because these are too small. Then we should ignore these chromosome arms
				if int(snpDivision[armInd]) < 1:
					emptyIndices.append(armInd)
			
		#	print len(self.allChromosomeArms)
			
			#Sometimes due to rounding not all SNPs are distributed. There should not be many so I append these to the last chromosome arm. 
			
			if len(self.allChromosomeArms) < self.snpNum:
				
				difference = self.snpNum - len(self.allChromosomeArms)
				
				#use the last arm which remains in memory
				self.allChromosomeArms += [self.chromosomeArms[armInd]]*difference
			
			
			
			
			#Remove the arms that do not have SNPs from the list of chromosome arms with SNPs. 
			tmpArms = np.array(self.chromosomeArms)
			self.chromosomeArms = list(np.delete(tmpArms, emptyIndices))
			
			#Also remove the arms without SNPs from the probabilities so tat these will not be chosen later in the simulations
			tmpArms = np.array(self.simulationProbabilities.armProbabilities)
			tmpArmsDeleted = np.delete(tmpArms, emptyIndices, axis=0)
			
			tupleFormattedArmProbabilities = []
			for arm in tmpArmsDeleted:
				tupleFormattedArmProbabilities.append(tuple(arm))
				
			self.simulationProbabilities.armProbabilities = tupleFormattedArmProbabilities	
			
		
		
		for arm in self.allChromosomeArms:
			
			
			
			self.chromosomes.append(arm) #or does this need to be without arm notation?
			
			armIndex = np.where(hg19CoordinateArray[:,0] == arm)[0]
			
			position = randint(hg19CoordinateArray[armIndex,1], (hg19CoordinateArray[armIndex,2]-1))
			self.positions.append(position)
			
	def getHg19Coordinates(self):
		
		hg19Coordinates = simulationSettings.files['segmentationFile']
		
		hg19CoordinateList = []
		with open(hg19Coordinates, 'rb') as f:
			for line in f:
				line = line.strip()
				
				splitLine = line.split("\t")
				chromosomeArm = splitLine[0]
				start = int(splitLine[1])
				end = int(splitLine[2])
				hg19CoordinateList.append([chromosomeArm, start, end])
				
		hg19CoordinateArray = np.array(hg19CoordinateList, dtype="object")
		
		return hg19CoordinateArray
		
	def distributeSnvMeasurementsAcrossGenome(self):
		
		#Given the start and end positions of chromosomes, randomly distribute SNVs across the genome.
		
		#1. Read the hg19 coordinates
		hg19CoordinateArray = self.getHg19Coordinates()
		
		#Make sure to remove arms that do not have any SNPs because otherwise we cannot link SNVs to chromosome arm loss/gains.
		hg19CoordinateArrayFiltered = []
		for arm in hg19CoordinateArray:
			
			if arm[0] in self.chromosomeArms:
				hg19CoordinateArrayFiltered.append(arm)
			
		hg19CoordinateArray = np.array(hg19CoordinateArrayFiltered)
		
		
		#Randomly choose positions for these SNVs.
		snvs = []
		#self.snvPositions = []
		for snv in range(0, self.snvNum):
			
			#Randomly choose an arm for this snv
			arm = randint(0, (hg19CoordinateArray.shape[0]-1))
			
			#Then determine a random position within this interval.
			position = randint(hg19CoordinateArray[arm,1], (hg19CoordinateArray[arm,2]-1))
			#self.snvPositions.append(position)
			snv = SomaticVariant(hg19CoordinateArray[arm,0], position, 0)
			snvs.append(snv)
		
		return snvs
	
	def obtainSomaticVariantIndices(self):
		
		print "number of positions: ", len(self.positions)
		print "arms: ", len(self.allChromosomeArms)
		
				
		offset = 0
		variantIndices = []
		for variant in self.snvs:
			
			position = variant.position
			chromosome = variant.chromosome

			#print "variant position: ", position

			for measurementPosition in range(0, len(self.positions)-1):
				
			#	print "measurement position: ", measurementPosition
			#	print "measurement position: ", self.positions[measurementPosition]
				
				
				if str(chromosome) == str(self.allChromosomeArms[measurementPosition]): #the chromosome needs to be the same
					print "match on chromosome", str(chromosome), str(self.allChromosomeArms[measurementPosition])
					#for all the variants within the SNP range		
					if int(position) > int(self.positions[measurementPosition]) and int(position) < int(self.positions[measurementPosition + 1]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
					#	print "adding offset", offset
						offset += 1
					#Situation where the somatic variant comes before the SNP measurements
					if str(chromosome) != str(self.allChromosomeArms[measurementPosition-1]) and int(position) <= int(self.positions[measurementPosition]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
					#	print "adding offset 2", offset
						offset += 1
					#Situation where the somatic variant comes after the SNP measurements
					if str(chromosome) != str(self.allChromosomeArms[measurementPosition+1]) and int(position) >= int(self.positions[measurementPosition]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
					#	print "adding offset 3", offset
						offset += 1
						
					#print "variant index: ", variantIndex
		
		print variantIndices
		#exit()
		
		return variantIndices
	
		
	def subclonalExpansion(self):

		pkl_file = open(simulationSettings.files['targetCloneInstance'], 'rb')
		targetClone = pickle.load(pkl_file)
		pkl_file.close()
		
		targetClone.segmentation = self.segmentation
		iteration = 1 #dummy, remove later
		#To do the horizontal shuffling, the same UUID must be provided for each run. So an sh script should run these scripts providing the UUIDs from the existing run. 
		#In case we wish to do horizontal permutations, check here if we need to generate the clones.
		#If not, skip this intial part and load the simulation data from the pkl file, shuffle the measurements, and then run TC.
		#Write the results to a different folder but with the same UUIDs
		if simulationSettings.runType['horizontalShuffle'] == True:
			newDir = simulationSettings.files['outputDir'] + self.uniqueID + '_horizontalShuffle'
		else:
			newDir = simulationSettings.files['outputDir'] + self.uniqueID
		
		os.makedirs(newDir)
		
		if simulationSettings.runType['horizontalShuffle'] == False:
			
			if simulationSettings.runType['mixedSamples'] == False:
				[samples, finalClones, realTree, savedMu] = self.generateSamples()
			else:
				[samples, finalClones, realTree, savedMu] = self.generateMixedSamples() #mixign subclones in a sample
				
			#obtain the C for every clone and merge it into a matrix
			
			cMatrix = np.empty([self.snpNum, len(finalClones)], dtype=float)
			for cloneInd in range(0, len(finalClones)):
				#make the C from the individual Cs
				cloneC = finalClones[cloneInd].C
				for cInd in range(0, len(cloneC)):
					cMatrix[cInd][cloneInd] = cloneC[cInd].c[1]
			
			#Also print A
					
			aMatrix = np.empty([self.snpNum, len(finalClones)], dtype=object)
			for cloneInd in range(0, len(finalClones)):
				cloneA = finalClones[cloneInd].A
				for aInd in range(0, len(cloneA)):
			
					aMatrix[aInd][cloneInd] = cloneA[aInd]
			
			#Check the somatic variants
			somVarMatrix = np.empty([len(finalClones[0].somaticVariants), len(finalClones)], dtype=float)
			for cloneInd in range(0, len(finalClones)):
				somVar = finalClones[cloneInd].somaticVariants
				for variant in range(0, len(somVar)):
					if somVar[variant].value == 1:
						somVarMatrix[variant][cloneInd] = 1
					else:
						somVarMatrix[variant][cloneInd] = 0
						
			
			#Keep the measurements to write to a file later
			measurementsMatrix = np.empty([self.snpNum, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].afMeasurements
				measurementsMatrix[:,cloneInd] = measurements
				
			lafMeasurementsMatrix = np.empty([self.snpNum, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].measurements.measurements
				lafMeasurementsMatrix[:,cloneInd] = measurements
			
			
			
			simulationDataHandler = SimulationDataHandler()
			simulationDataHandler.setSimulationData(cMatrix, aMatrix, samples, measurementsMatrix, lafMeasurementsMatrix, somVarMatrix, realTree, self.chromosomes, self.positions)
			
			with open(newDir + '/simulationData.pkl', 'wb') as handle:
				pickle.dump(simulationDataHandler, handle, protocol=pickle.HIGHEST_PROTOCOL)
		else:
			
			#Load the data from the pkl file
			[simulationData, savedMu] = self.parseSimulationData(self.uniqueID)
			
			cMatrix = simulationData.cMatrix
			aMatrix = simulationData.aMatrix
			realTree = simulationData.realTree
			samples = simulationData.samples
			measurementsMatrix = simulationData.afMatrix
			lafMeasurementsMatrix = simulationData.lafMatrix
			
			#The AF need to be shuffled randomly within each sample
			
			for sampleInd in range(0, len(samples)):
				sample = samples[sampleInd]
				
				afMeasurements = sample.afMeasurements
				lafMeasurements = sample.measurements.measurements
				
				#Determine the random positions that the new measurements will have
				newOrder = random.sample(range(0, len(afMeasurements)), len(afMeasurements))
				
				
				shuffledAFMeasurements = []
				shuffledLAFMeasurements = []
				for measurementInd in newOrder:
					
					shuffledAFMeasurements.append(afMeasurements[measurementInd])
					shuffledLAFMeasurements.append(lafMeasurements[measurementInd])
			
				
				samples[sampleInd].afMeasurements = shuffledAFMeasurements
				samples[sampleInd].measurements = self.generateLAFObject(shuffledLAFMeasurements)
			
			somVarMatrix = simulationData.snvMatrix
			
			measurementsMatrix = np.empty([self.snpNum, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].afMeasurements
				measurementsMatrix[:,cloneInd] = measurements
				
			lafMeasurementsMatrix = np.empty([self.snpNum, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].measurements.measurements
				lafMeasurementsMatrix[:,cloneInd] = measurements
		
		
			
		#Run TC	
		[eCMatrix, eAMatrix, eSamples, trees, iterationMu, iterationMessages] = targetClone.run(samples)
		
		
		#1. Implement different scoring for the output
		simulationErrorHandler = SimulationErrorHandler()
		
		allCScores = []
		allAScores = []
		allMuScores = []
		allTreeScores = []
		allAmbiguityScores = []
		allAmbiguityCorrectedScores = []
		
		
		eCMatrixFloat = eCMatrix.astype(float)
		cMatrixFloat = cMatrix.astype(float)
		#cScore = self.computeCRMSE(cMatrixFloat, eCMatrixFloat)
		cScore = simulationErrorHandler.computeCError(cMatrixFloat, eCMatrixFloat)
		
		allCScores.append(cScore / cMatrixFloat.size)
		
		#aScore = self.computeARMSE(aMatrix, eAMatrix)
		aData = simulationErrorHandler.computeAError(aMatrix, eAMatrix)
		aScore = aData[0] / float(aMatrix.size)
		
		allAScores.append(aScore)
		eAStringMatrix = aData[2]
		aStringMatrix = aData[1]
		
		muData = simulationErrorHandler.computeMuError(eSamples, savedMu)
		muScore = muData[0]
		allMuScores.append(muScore)
		allRealMuT = muData[1]
		allEMuT = muData[2]
		
		treeScore = simulationErrorHandler.computeTreeError(trees, realTree)
		allTreeScores.append(treeScore)
		bestTree = trees[len(trees)-1]
		#for each iteration we can save the individual scores. We use these later to make some plots of the current performance.
		
		#Also save the simulated data to files, we can then look back at these later if we need to have more information.
	
		#Somatic variants matrix
		np.savetxt(newDir + '/SomVar_' + str(iteration) + '.txt', somVarMatrix, fmt='%i', delimiter='\t')
					

		#Save measurements
		np.savetxt(newDir + '/AFMeasurements_' + str(iteration) + '.txt', measurementsMatrix, fmt='%f', delimiter='\t')
		np.savetxt(newDir + '/LAFMeasurements_' + str(iteration) + '.txt', lafMeasurementsMatrix, fmt='%f', delimiter='\t')
		
		#Save C
		np.savetxt(newDir + '/RealC_' + str(iteration) + '.txt', cMatrix.astype(int), fmt='%i', delimiter='\t')
		np.savetxt(newDir + '/EstimatedC_' + str(iteration) + '.txt', eCMatrix.astype(int), fmt='%i', delimiter='\t')
		
		#Save A
		np.savetxt(newDir + '/RealA_' + str(iteration) + '.txt', aStringMatrix, fmt='%s', delimiter='\t')
		np.savetxt(newDir + '/EstimatedA_' + str(iteration) + '.txt', eAStringMatrix, fmt='%s', delimiter='\t')
		
		#Save Mu
		np.savetxt(newDir + '/RealMu_' + str(iteration) + '.txt', allRealMuT, fmt='%f', delimiter='\t')
		np.savetxt(newDir + '/EstimatedMu_' + str(iteration) + '.txt', allEMuT, fmt='%f', delimiter='\t')
		
		#Save the trees
		f = open(newDir + '/RealTrees_' + str(iteration) + '.txt', 'w')
		f.write(str(realTree.getGraph()))  # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/EstimatedTrees_' + str(iteration) + '.txt', 'w')
		f.write(str(bestTree.getGraph()))  # python will convert \n to os.linesep
		f.close()
	
	
		[ambiguityScore, totalScore] = self.computeAmbiguityScore(aMatrix, eAMatrix, allRealMuT, allEMuT)
		allAmbiguityScores.append(ambiguityScore)
		allAmbiguityCorrectedScores.append(totalScore)
		print "ambiguity error: ", ambiguityScore
		
		
		
		print "all C errors: ", allCScores
		
		
		print "all A errors: ", allAScores
		
		
		print "all mu errors: ", allMuScores
		
		
		print "all Tree errors: ", allTreeScores
		
		f = open(newDir + '/cError.txt', 'w')
		
		f.write(str(allCScores[0]))  # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/aError.txt', 'w')
		f.write(str(allAScores[0]))  # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/muError.txt', 'w')
		f.write(str(allMuScores[0]))  # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/treeError.txt', 'w')
		f.write(str(allTreeScores[0])) # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/ambiguityError.txt', 'w')
		f.write(str(allAmbiguityScores[0])) # python will convert \n to os.linesep
		f.close()
		f = open(newDir + '/ambiguityCorrectedError.txt', 'w')
		f.write(str(allAmbiguityCorrectedScores[0])) # python will convert \n to os.linesep
		f.close()
		
	def readDataFromFile(self, file):
		text_file = open(file, "r")
		lines = text_file.read()
		floatLines = []
		
		for line in lines.split("\n"):
			if line != "":
				floatLines.append(float(line))
		
		text_file.close()
		
		return floatLines
	
	def parseSimulationData(self, uniqueIdentifier):
		#Load the file by unique identifier, and store the object
		#Get the actual location from the settings file!
		pkl_file = open(simulationSettings.files['outputDir'] + uniqueIdentifier + '/simulationData.pkl', 'rb')
		simulationData = pickle.load(pkl_file)
		pkl_file.close()
		
		
		#We will have to obtain the real mu from the file.
		muFile = simulationSettings.files['outputDir'] + uniqueIdentifier + '/RealMu_1.txt'
		muData = self.readDataFromFile(muFile)
		
		
		savedMu = []
		for muInd in range(0, len(muData)):
			mu = muData[muInd]
			
			muObject = Mu(1 - float(mu)) #the mu in the file is the tumor mu, we need to provide the normal cell mu here!
			savedMu.append(muObject)
		
		
		return [simulationData, savedMu]
		
	def generateSamples(self):
		
		#Starting from a healthy cell, start making subclones	
			
		healthySubclone = Subclone()
		healthySubclone.C = []
		
		for snpInd in range(0, self.snpNum):
			healthySubclone.C.append(C([2,2]))
		

		healthySubclone.name = 'Healthy'
		
		healthyAlleleObjects = []
		prevArm = self.allChromosomeArms[0]
		alleleA = 'A' + str(uuid.uuid4())
		alleleB = 'B' + str(uuid.uuid4())
		
		for newAlleleInd in range(0, self.snpNum):
			
			if self.allChromosomeArms[newAlleleInd] != prevArm:
				alleleA = 'A' + str(uuid.uuid4())
				alleleB = 'B' + str(uuid.uuid4())
				
			newAlleles = Alleles(1,1)
			newAlleles.alleleIdentifiers = [alleleA, alleleB] #we add alleles A and B
			healthyAlleleObjects.append(newAlleles)
			prevArm = self.allChromosomeArms[newAlleleInd]
		
		healthySubclone.A = healthyAlleleObjects
		
		healthySubclone.somaticVariants = self.snvs
		
		#Make a sample object for this subclone
		healthySample = Sample(None, None)
		healthySample.C = healthySubclone.C
		healthySample.A = healthySubclone.A
		healthySample.Mu = Mu(100)
		#obtain the chromosome, start and end information from the other samples
		measurements = self.generateMeasurements(healthySubclone, 0) #the mu of the tumor is here 0 because we have a healthy sample
		healthySample.afMeasurements = measurements[0]
		healthySample.measurements = measurements[1]

		healthySample.somaticVariants = [0]*self.snvNum #Why is this not the same as the object-based way?
		healthySample.somaticVariantsInd = self.snvPositions
		healthySample.setParent(None)
		healthySample.name = healthySubclone.name
		
		#Do we still need these values? Add them anyway for now
		eventDistances = EventDistances(self.kmin, self.kmax)
		bestCMuHealthy = CMuCombination(C([2,2]), Mu(100), eventDistances)
		
		healthySample.bestCMu = [bestCMuHealthy]*len(measurements[1].measurements)
		healthySample.originalCMu = healthySample.bestCMu

	
		#Then do the subclonal expansion starting from this healthy subclone
		clones = self.evolve(healthySubclone, 1, 1, self.simulationProbabilities, [])
		finalClones = [healthySubclone]
		cloneInd = 0
		savedMu = [healthySample.Mu]
		edgeList = []
		samples = [healthySample]
		sampleNames = []
		for clone in clones: #make a new sample for each clone
			#if cloneInd in sampleInd:
			finalClones.append(clone)
			newSample = Sample(None, None)
			#We cannot set C, A and mu because we do not know it yet
			
			somVarArray = []
			for variant in clone.somaticVariants:
				if variant.value == 1:
					somVarArray.append(1)
				else:
					somVarArray.append(0)
			
			
			newSample.somaticVariants = somVarArray
			newSample.somaticVariantsInd = self.snvPositions
			
			randomMuT = random.sample(range(0, 100), 1)[0] #mu of the tumor component
			
			newSampleMu = Mu(self.muN)
			
			savedMu.append(newSampleMu)
			#Calculate these measurements from the counts
			measurements = self.generateMeasurements(clone, 100-self.muN)
			newSample.afMeasurements = measurements[0]
			newSample.measurements = measurements[1]
			
			#newSample.setParent(healthySample) #we do not know the parent, this will first be the healthy sample for everyone
			#Do not update the sample information just yet, as we wish to still define the real underlying tree first
			newSample.name = clone.name

			parentName = clone.parent.name
			newEdge = (0,parentName,newSample.name) #use a weight of 0 by default, this is alright for comparison, we do not know the actual weight
			edgeList.append(newEdge)
			
			samples.append(newSample)
			sampleNames.append(newSample.name)
			cloneInd += 1
		
		vertices = []
		for edge in edgeList:
			parent = edge[1]
			child = edge[2]
			if parent not in vertices:
				vertices.append(parent)
			if child not in vertices:
				vertices.append(child)
		
		realTree = Graph(vertices, set(edgeList), edgeList)
		
		
		if simulationSettings.runType['randomMeasurements'] == True: #assign random SNV measurements to the samples
			varCount = len(samples[0].somaticVariants)
			
			for sample in samples:
				
				#Make an array with random SNVs
				randomSNVs = []
				for newVar in range(0, varCount):
					randomValue = randint(0,1)
					randomSNVs.append(randomValue)
				sample.somaticVariants = randomSNVs
		
		return samples, finalClones, realTree, savedMu
	
	#Function to generate samples but then when multiple subclones are mixed within a sample. 
	def generateMixedSamples(self):
		#Starting from a healthy cell, start making subclones	
			
		healthySubclone = Subclone()
		healthySubclone.C = []
		
		for snpInd in range(0, self.snpNum):
			healthySubclone.C.append(C([2,2]))
		

		healthySubclone.name = 'Healthy'
		
		healthyAlleleObjects = []
		prevArm = self.allChromosomeArms[0]
		alleleA = 'A' + str(uuid.uuid4())
		alleleB = 'B' + str(uuid.uuid4())
		
		for newAlleleInd in range(0, self.snpNum):
			
			if self.allChromosomeArms[newAlleleInd] != prevArm:
				alleleA = 'A' + str(uuid.uuid4())
				alleleB = 'B' + str(uuid.uuid4())
				
			newAlleles = Alleles(1,1)
			newAlleles.alleleIdentifiers = [alleleA, alleleB] #we add alleles A and B
			healthyAlleleObjects.append(newAlleles)
			prevArm = self.allChromosomeArms[newAlleleInd]
		
		healthySubclone.A = healthyAlleleObjects
		
		healthySubclone.somaticVariants = self.snvs
		
		#Make a sample object for this subclone
		healthySample = Sample(None, None)
		healthySample.C = healthySubclone.C
		healthySample.A = healthySubclone.A
		healthySample.Mu = Mu(100)
		#obtain the chromosome, start and end information from the other samples
		measurements = self.generateMeasurements(healthySubclone, 0) #the mu of the tumor is here 0 because we have a healthy sample
		healthySample.afMeasurements = measurements[0]
		healthySample.measurements = measurements[1]

		healthySample.somaticVariants = [0]*self.snvNum #Why is this not the same as the object-based way?
		healthySample.somaticVariantsInd = self.snvPositions
		healthySample.setParent(None)
		healthySample.name = healthySubclone.name
		
		#Do we still need these values? Add them anyway for now
		eventDistances = EventDistances(self.kmin, self.kmax)
		bestCMuHealthy = CMuCombination(C([2,2]), Mu(100), eventDistances)
		
		healthySample.bestCMu = [bestCMuHealthy]*len(measurements[1].measurements)
		healthySample.originalCMu = healthySample.bestCMu

	
		#Then do the subclonal expansion starting from this healthy subclone
		clones = self.evolve(healthySubclone, 1, 1, self.simulationProbabilities, [])
		finalClones = [healthySubclone]
		cloneInd = 0
		savedMu = [healthySample.Mu]
		edgeList = []
		samples = [healthySample]
		sampleNames = []
		
		numberOfMixedClones = simulationSettings.general['numberOfMixedClones'] #Number of subclones to mix within a sample
		maximumMinorCloneFrequency = simulationSettings.general['maximumMinorCloneFrequency']
		minimumMinorCloneFrequency = simulationSettings.general['minimumMinorCloneFrequency']
		
		print "total number of clones: ", len(clones)
		
		for sample in range(0, len(clones)):
			
			#Obtain the main clone information
			clone = clones[sample]
			print "processing clone: ", clone.name
			finalClones.append(clone)
			
			
			possibleSamplesToMix = range(0,sample) + range(sample+1,len(clones)) #the samples to select can never be the sample itself
			print "can mix with: ", possibleSamplesToMix
			print "number of mixed clones: ", numberOfMixedClones
			
			#If our tree does not contain enough subclones to mix with, we mix with the maximum allowed.
			numberOfClonesToMix = numberOfMixedClones
			if numberOfClonesToMix >= len(clones):
				numberOfClonesToMix = len(clones)-1
			
			subclonesToMix = random.sample(possibleSamplesToMix, numberOfClonesToMix) #Select x other samples to mix the current sample with
			
			print "mixing with: ", subclonesToMix
			
			selectedClones = [clone] #make sure that we use this sample as major clone
			for subcloneInd in subclonesToMix:
				selectedClones.append(clones[subcloneInd])
			
			
			#We wish to have a maximum value of 10 for the minor subclones. This function below makes vectors with a sum of 10. Thus, we can sample 1 more value than we need
			#and then randomly select n-1 from this vector to get two values between 1 and 10. Make sure that this value that we sample is not 0!
			
			print maximumMinorCloneFrequency,numberOfMixedClones+1,maximumMinorCloneFrequency
			
			#Keep sampling if the value is smaller than the minimum clone fraction (maybe a bad solution but I don't have a better one right now)
			totalMinorFraction = 0
			
			while totalMinorFraction < minimumMinorCloneFrequency:
			
				#sum to equal, N to sample, maxval, minval
				minorCloneFractionsFull = randvec(maximumMinorCloneFrequency,numberOfMixedClones+1,maximumMinorCloneFrequency,1) #Sample the minor clone fractions
				
				print "minor fractions: ", minorCloneFractionsFull
				minorCloneFractions = random.sample(minorCloneFractionsFull, numberOfMixedClones)
				
				totalMinorFraction = sum(minorCloneFractions)
				print "minor fractions filtered: ", minorCloneFractions
				print "total minor fraction: ", totalMinorFraction
				
			
			#Do +1 to make sure that the major clone fraction is also not the same as the minor clone fraction
			
			normalFraction = 10 #test to see influence of the normal fraction
			
			majorFraction = [100 - normalFraction - totalMinorFraction]
			print "major fraction: ", majorFraction
			print "minor fraction: ", totalMinorFraction
			print minorCloneFractions
			#Combine the minor fractions with the major fractions
			mixingFractions = majorFraction + minorCloneFractions
			
			print "mixing fractions: ", mixingFractions
			
			tumorFraction = sum(mixingFractions)
			normalFraction = 100 - tumorFraction
			print "total tumor fraction: ", tumorFraction
			print "normal fraction: ", normalFraction
			
			allMixingFractions = dict()
			allMixingFractions[sample] = mixingFractions
			
			#Generate AF measurements and somatic SNV measurements
			
			sampleAFMeasurements = self.generateMixedSampleAFMeasurements(selectedClones, mixingFractions)
			sampleSNVMeasurements = self.generateMixedSampleSNVMeasurements(selectedClones, mixingFractions)
			#Do not mix SNVs, now use only the SNVs of the original clone
			# cloneSnvMeasurements= selectedClones[0].somaticVariants
			# sampleSNVMeasurements = []
			# for snv in cloneSnvMeasurements:
			# 	sampleSNVMeasurements.append(snv.value)
			# #sampleSNVMeasurements = selectedClones[0].somaticVariants
			# print sampleSNVMeasurements
			
			#Generate LAF measurements from AF measurements
			sampleLAFMeasurements = self.generateLAFObject(sampleAFMeasurements)
			
			
			#Create a sample object
			newSample = Sample(None, None)
			#We cannot set C, A and mu because we do not know it yet
	
			newSample.somaticVariants = sampleSNVMeasurements
			newSample.somaticVariantsInd = self.snvPositions
			
			newSampleMu = Mu(100 - sum(mixingFractions)) #Set the fraction of the normal component
			
			savedMu.append(newSampleMu)
			#Calculate these measurements from the counts
			newSample.afMeasurements = sampleAFMeasurements
			newSample.measurements = sampleLAFMeasurements
			
			
			
			#The relations between the samples are no longer the same, since clones have been mixed. 
			newSample.name = clone.name
		
			parentName = clone.parent.name
			newEdge = (0,parentName,newSample.name) #use a weight of 0 by default, this is alright for comparison, we do not know the actual weight
			edgeList.append(newEdge)
			
			samples.append(newSample)
			sampleNames.append(newSample.name)

	
		vertices = []
		for edge in edgeList:
			parent = edge[1]
			child = edge[2]
			if parent not in vertices:
				vertices.append(parent)
			if child not in vertices:
				vertices.append(child)
		
		realTree = Graph(vertices, set(edgeList), edgeList)
		
		return samples, finalClones, realTree, savedMu

	
	def generateMixedSampleSNVMeasurements(self, clones, mixingFractions):

		#We should compute the VAF that the SNV would have.
		
		#I think that the positions that we use are before/after. This is thus only an estimate!
		
		#Per clone, we count the total A and B fraction for that SNV.
		#A fraction is the total number of alleles at that locus - the number of alleles with an SNV
		#B fraction is the number of alleles with an SNV
		#This should be multiplied with the mixing fraction of the clone
		#Normal cells have an A and B fraction of 0.
		
		snvMeasurements = []
		#Obtain the A and B count at each position in the AF measurements
		
		#For position, for clone:
		#Obtain A and B Count
		#Multiply with the right mixing fraction
		#store as AF measurments
		
		variants = len(clones[0].somaticVariants)
		
		for variant in range(0, variants):
			totalRefFraction = 0
			totalVarFraction = 0
			
			#Also mix in the normal component (will always be 0)
			normalRefCount = 2 
			normalVarCount = 0
			normalFraction = 100 - sum(mixingFractions)
			
			totalRefFraction += normalRefCount * normalFraction
			totalVarFraction += normalVarCount * normalFraction
			
			for cloneInd in range(0, len(clones)):
				clone = clones[cloneInd]
				
				
				
				#The allele is the same across the chromosome arm, so just look anywhere at that chromosome arm
				somVarChromosome = clone.somaticVariants[variant].chromosome
				
				chromosomeIndex = self.allChromosomeArms.index(somVarChromosome) #get the first index of a snp position on this chromosome
				
				#Obtain the allele as a string to count how many are A and B (is probably just the first and second entry of array, could be easier than this)
				allele = clone.A[chromosomeIndex].getAllelesAsString()
				#print "alleles of clone: ", allele
				#Obtain how many A and B alleles are present
				AOccurrences = [m.start() for m in re.finditer('A', allele)]
				ACount = len(AOccurrences)
				
				BOccurrences = [m.start() for m in re.finditer('B', allele)]
				BCount = len(BOccurrences)
				
				#Compute the total alleles
				totalAlleles = ACount + BCount
				
				#Obtain the number of alleles the SNV is present on
				snvAlleles = len(clone.somaticVariants[variant].alleles)
				
				# print "current clone: ", cloneInd
				# print "mixing fractions: ", mixingFractions
				# print "total alleles: ", totalAlleles
				# print "snv alleles: ", snvAlleles
				
				refAlleles = totalAlleles - snvAlleles
				
				#Mix the allele counts with mixing fractions to compute the fraction of A and B alleles in this sample given the mixing fractions
				refFraction = refAlleles * mixingFractions[cloneInd]
				varFraction = snvAlleles * mixingFractions[cloneInd]
				
				#Sum the fractions to eventually divide the fraction of the B allele (or A) by the total.
				totalRefFraction += refFraction
				totalVarFraction += varFraction
			
			totalAlleleCount = totalRefFraction + totalVarFraction

			#Compute the VAF
			vaf = totalVarFraction / float(totalAlleleCount)
			
			#Do a filter based on the settings to filter out SNVs from contamination. 0.04 may not be enough!
			#Then print to file and see how the SNV patterns look. Does it improve tree reconstruction?
			if vaf < simulationSettings.general['minimumSNVFrequency']:
				#print "true"
				vaf = 0 #we consider these SNVs as not present when their frequency is too low. 
			else:
				vaf = 1
			snvMeasurements.append(vaf)
	
		return snvMeasurements
		
	def generateMixedSampleAFMeasurements(self, clones, mixingFractions):
		afMeasurements = []
		#Obtain the A and B count at each position in the AF measurements
		
		#For position, for clone:
		#Obtain A and B Count
		#Multiply with the right mixing fraction
		#store as AF measurments
		
		positions = len(clones[0].A)
		
		for position in range(0, positions):
			totalAFraction = 0
			totalBFraction = 0
			#totalAlleleCount = 0
			
			#Also mix in the normal component
			normalACount = 1
			normalBCount = 1
			normalFraction = 100 - sum(mixingFractions)
			
			totalAFraction += normalACount * normalFraction
			totalBFraction += normalBCount * normalFraction
			
			for cloneInd in range(0, len(clones)):
				clone = clones[cloneInd]
				#Obtain the allele as a string to count how many are A and B (is probably just the first and second entry of array, could be easier than this)
				allele = clone.A[position].getAllelesAsString()
				#print allele
				#Obtain how many A and B alleles are present
				AOccurrences = [m.start() for m in re.finditer('A', allele)]
				ACount = len(AOccurrences)
				
				BOccurrences = [m.start() for m in re.finditer('B', allele)]
				BCount = len(BOccurrences)
				
				#Mix the allele counts with mixing fractions to compute the fraction of A and B alleles in this sample given the mixing fractions
				aFraction = ACount * mixingFractions[cloneInd]
				bFraction = BCount * mixingFractions[cloneInd]
				
				#Sum the fractions to eventually divide the fraction of the B allele (or A) by the total.
				totalAFraction += aFraction
				totalBFraction += bFraction
			
			totalAlleleCount = totalAFraction + totalBFraction

			#Compute the AF
			af = totalBFraction / float(totalAlleleCount)
			#Add noise to the AF
			lower, upper = 0, 1 #AF truncated ones
			mu, sigma = af, simulationSettings.general['noiseLevel']
			noisedAf = af
			if simulationSettings.general['noiseLevel'] > 0:
				X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
				noisedAf = X.rvs(1)[0] #take a sample
			
			
				
			
			afMeasurements.append(noisedAf)
	
		return afMeasurements
	
	


	def computeAmbiguityScore(self, aMatrix, eAMatrix, realMuT, eMuT):
		#define some random matrices to test
		ambiguityError = 0
		totalError = 0
		for row in range(0, aMatrix.shape[0]):
			ambiguityFound = False #check if we find an ambiguity somewhere on this line. If we do, the row error should be 0 (we do not count errors)
			rowError = 0
			for col in range(0, aMatrix.shape[1]):
				currentMu = realMuT[col] #mu of the current subclone
				currentEstimatedMu = eMuT[col]
				realA = aMatrix[row][col]
				#for the real A, compute the expected LAF.
				aCountT = realA.ACount * currentMu
				bCountT = realA.BCount * currentMu
				
				aCountN = 1 * (1 - currentMu)
				bCountN = 1 * (1 - currentMu)
				
				totalACount = aCountT + aCountN
				totalBCount = bCountT + bCountN
				#the LAF is the minimum of the a and b count / the sum of both
				
				countSum = totalACount + totalBCount
				minCount = min([totalACount] + [totalBCount])
				
				realLaf = minCount / countSum
				
				#repeat this for the estimated a
				estimatedA = eAMatrix[row][col]
				#for the real A, compute the expected LAF.
				aCountT = estimatedA.ACount * currentEstimatedMu
				bCountT = estimatedA.BCount * currentEstimatedMu
				
				aCountN = 1 * (1 - currentEstimatedMu)
				bCountN = 1 * (1 - currentEstimatedMu)
				
				totalACount = aCountT + aCountN
				totalBCount = bCountT + bCountN
				#the LAF is the minimum of the a and b count / the sum of both
				
				countSum = totalACount + totalBCount
				minCount = min([totalACount] + [totalBCount])
				
				estimatedLaf = minCount / countSum
				
				if realA.ACount != estimatedA.ACount or realA.BCount != estimatedA.BCount: #if the alleles are different and the LAF is also different, error.
					if ambiguityFound == False:
						rowError += 1
				
				#if the laf are the same, but the A is different, add to the ambiguity score.
				#otherwise add to the error.
				if realLaf == estimatedLaf:
					if realA.ACount != estimatedA.ACount or realA.BCount != estimatedA.BCount: #the alleles are different if at least the A or B score is different. One of them can be the same
						ambiguityError += 1
						ambiguityFound = True
						rowError = 0
						#if we find an ambiguity somewhere on this row, the rowError should be 0 and not increase. 
			totalError += rowError 		
		return [(ambiguityError / float(aMatrix.size)), (totalError / float(aMatrix.size))]		
	
	def evolve(self, subclone, cycle, round, simulationProbabilities, clones):
		print "current cycle: ", cycle
		
		if cycle == self.cellCycles or cycle > self.cellCycles: 
			return clones
		
		#Otherwise, evolve this subclone into a new one and evolve it until we are done
		
		#newSubclone = subclone.expandSubclone(simulationProbabilities, self)
		#we have two branches: we start with the malignant one and let it develop into two subclones.
		#The first is the new subclone, the second is the original.
		#the first goes back into the function and is further evolved
		#then the second also goes back into the function.
		
		#Make a new subclone, expand it from the previous and save the result. 
		newSubclone = subclone.expandSubclone(simulationProbabilities, self)
		
		#print "original clone: ", subclone
		#print "new subclone: ", newSubclone
		if newSubclone == False: #do not expand the clones further if these are not viable.
			print "unviable"
			
			passed = False
			currentCycle = 0
			if subclone.name == 'GCNIS': #if the parent is GCNIS, we assume that we are unlimited in our sources and we can keep searching until we find a viable subclone. 
				print "searching for viable clones"
				while passed == False:
					print "current cycle: ", currentCycle
					newSubclone = subclone.expandSubclone(simulationProbabilities, self)
					currentCycle += 1
					if newSubclone != False:
						print "viable at cycle ", currentCycle
						passed = True
			else:
				return clones
			
		else:
			#clones.append(subclone)
			clones.append(newSubclone) #only append the new subclone when it is viable
		
		clones = self.evolve(subclone, cycle + 1, round, simulationProbabilities, clones) #left branch
		#leftMostInd += 1
		clones = self.evolve(newSubclone, cycle + 1, round, simulationProbabilities, clones) #right branch
		
		return clones

	def generateLAFObject(self, measurements):
		
		return LAF(measurements, self.chromosomes, self.positions, self.positions)
		
	def generateMeasurements(self, sample, randomMuT):
		
		if simulationSettings.runType['randomMeasurements'] == True: #in this case sample random LAF
		
			#Generate random LAF and AF measurements
			
			AF = []
			LAF = []
			
			for alleles in sample.A:
				randomAF = randint(0,1000)
				randomAF = randomAF / float(1000)
				
				randomLAF = randomAF
				if randomAF > 0.5:
					randomLAF = round(1 - randomAF, 3)
				AF.append(randomAF)
				LAF.append(randomLAF)

		else:
			
			#The LAF is easily computed for every position from the A matrix
			randomMuT = randomMuT / float(100)
			#The values here need to be weighed with a random mu, let's generate a random one here and see how it works out.
			#we also need to add the correct number of alleles for healthy cells to generate the correct AFs
			
			#Generate the measurement values from a normal distribution with added noise (first do this without noise)
			
			AF = []
			LAF = []
			#currentNoiseLevel = self.noiseLevels[0]
			
			#currentNoiseLevel = 0
			currentNoiseLevel = simulationSettings.general['noiseLevel']
			#For each position that we can add a measurement to
			#Check how many chromosome names are there, add this many measurements (later we add random noise at this step)
			#print "arms: ", self.possibleArms
			
			measurementPosition = 0
			#print "carms: ", len(self.chromosomeArms)
			#print "number of measurements: ", len(sample.A)
			for alleles in sample.A:
				#print measurementPosition
				#What is the number of measurements on this chromosome arm?
				#Add the same number of measurements.
	
				currentArm = self.allChromosomeArms[measurementPosition] #at this point sample does not know where the measurements are as this is part of the measurement object
				armIndices = [i for i, x in enumerate(self.allChromosomeArms) if x == currentArm]
				numberOfMeasurements = len(armIndices)
	
				aCountT = alleles.ACount * randomMuT
				bCountT = alleles.BCount * randomMuT
				
				aCountN = 1 * (1 - randomMuT)
				bCountN = 1 * (1 - randomMuT)
				
				totalACount = aCountT + aCountN
				totalBCount = bCountT + bCountN
				#the LAF is the minimum of the a and b count / the sum of both
				countSum = totalACount + totalBCount
				minCount = min([totalACount] + [totalBCount])
				if countSum == 0:
					newAf = 0
				else:
					#Add noise to the measurements
					newAf = totalBCount / float(countSum) #we assume that the b count is the one that is the variant allele. We do not know!
				
				#take a sample from a normal distribution where the mu is the AF position (or should we sample from the LAF?)
				#noisedAf = np.random.normal(newAf, currentNoiseLevel, 1)[0]
				
				if currentNoiseLevel > 0:
					
					lower, upper = 0, 1 #AF truncated ones
					mu, sigma = newAf, currentNoiseLevel
					X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
					noisedAf = X.rvs(1)[0] #take a sample
					#noisedAf = newAf
					#AF.append(noisedAf)
					
				else:
					noisedAf = newAf
					
				AF.append(noisedAf)
				if noisedAf > 0.5:
					LAF.append(1-noisedAf)
				else:
					LAF.append(noisedAf)
					
				#for measurementNum in range(0, numberOfMeasurements): #here we should later add noise. 
				#LAF.append(minCount / float(countSum))
				
				measurementPosition += 1
		
		return [AF, self.generateLAFObject(LAF)]
		
		
	

class SimulationProbabilities:
	
	somaticVariants = None #this can be a list with objects
	armProbabilities = None #make this a list of tuples to preserve sorting. 
		
	#read the file with probabilities
	
	def readFullProbabilityFile(self, file):
		self.armProbabilities = []
		self.somaticVariants = []
		#Read the file, store the somatic variant positions and arm probabilities separately
		varCounter = 0
		with open(file, "r") as inFile:
			array = []
			lineNum = 0
			previousLine = ''
			section = ''
			for line in inFile:
				line = line.strip('\r\n') #remove any possible newlines
				
				#Check if this is a header and which section we are in
				if re.match("^\#\#.*", line):
					if section == '':
						section = 'somVar'
					else:
						section = 'chrArms'
					previousLine = 'header'
					continue
				
				#if we end up here, we are reading the column names
				if previousLine == 'header':
					previousLine = ''
					continue
				
				#otherwise, this is a line with interesting data.
				splitLine = re.split("\t", line)
				if section == 'somVar':
					somVar = SomaticVariant(splitLine[0], splitLine[1], 0)
					somVar.ind = varCounter
					self.somaticVariants.append(somVar)
					varCounter += 1
				else:
					#self.armProbabilities[splitLine[0]] = (splitLine[1], splitLine[2].replace(",", "."), splitLine[3].replace(",", "."))
					self.armProbabilities.append((splitLine[0],splitLine[1], splitLine[2].replace(",", "."), splitLine[3].replace(",", ".")))
		
		
	#store the probabilities separately
	#if we wish to introduce a new event, we ask this object for a certain event and it will return it to us.
	#for the somatic variants, we need to remember which once have already been used before.
	#the
	
	def readProbabilities(self, file):
		#Read the probabilities from the file, these are in columns 2 and 3
		
		with open(file) as tsv:
			colNum = 0
			probabilities = None
			for column in zip(*[line for line in csv.reader(tsv, dialect="excel-tab")]):
				if colNum == 0:
					probabilities = np.zeros((len(column),2))
				if colNum == 2:
					column = [float(i.replace(",", ".")) for i in list(column)]
					
					probabilities[:,0] = column
					#store the values in an array
				if colNum == 3:
					column = [float(i.replace(",", ".")) for i in list(column)]
					probabilities[:,1] = column
				colNum += 1
	
		return probabilities
	
	def minMaxNormalizedProbabilities(self, probabilities):
		
		x = 0
		y = 1
		minProbabilities = np.min(probabilities)
		maxProbabilities = np.max(probabilities)
		probabilityRange = maxProbabilities - minProbabilities
		normalizedProbabilities = (probabilities - minProbabilities) / float(probabilityRange)
		
		probabilityHeightRange = y - x
		normalizedProbabilities = (normalizedProbabilities*probabilityHeightRange) + x
		
		return normalizedProbabilities
	
	def normalizeProbabilitiesPerArm(self, probabilities):
		
		normalizedProbabilitiesPerArm = np.zeros(probabilities.shape)
		
		for rowInd in range(0, probabilities.shape[0]):
			#obtain the values of this row
			rowValues = probabilities[rowInd,:]
			#compute the sum of both values
			summedValues = np.sum(rowValues)
			#print summedValues
			normalizedRow = rowValues / summedValues
			#print normalizedRow
			normalizedProbabilitiesPerArm[rowInd,:] = normalizedRow
		return normalizedProbabilitiesPerArm	

	def convertProbabilities(self, file, outfile):
		#The original probabilities are not useful for us:
		#Read the original file
		probabilities = self.readProbabilities(file)
	
		#Do a min-max normalization across all loss and gain proobabilities.
		normalizedProbabilities = self.minMaxNormalizedProbabilities(probabilities)
		#Then normalize the data per chromosome arm.
		normalizedProbabilitiesPerArm = self.normalizeProbabilitiesPerArm(normalizedProbabilities)
		
		#Write these data to a new file
		np.savetxt(outfile, normalizedProbabilities, fmt='%f', delimiter='\t')

	
simulator = Simulator()



	