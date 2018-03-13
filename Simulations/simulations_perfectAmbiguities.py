#Classes for generating simulation data

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
from permutations import Permute


class Simulator:
	
	#Define the number of positions that can be affected by SVs, losses/gains and amplifications in the evolution process
	#The exact definitions of these events are defined in an excel file. 
	numberOfSVs = 36
	numberOfArmChanges = 20
	numberOfAllelicAmps = 4
	
	#Define the number of changes that can happen during each simulated cell cycle
	malignantWcLosses = 10
	malignantArmLosses = 10
	malignant12pCopies = 6
	malignantSVGains = 20
	
	cycleSVGains = 2
	cycleArmGains = 3
	cycleArmLosses = 8
	
	kmin = 1
	kmax = 6
	cellCycles = 4
	
	referenceFile = 'reference.txt' #In reference.txt we have stored the selected positions from T3209 that are not 0 in the reference sample. We have the chromosome numbers and positions. From the
	#measured AF values we estimate a SD that we use to generate a noise distribution from which we sample AF measurements. We generate new measurements at these specified positions. 
	segmentationFile = 'pq_segmentation.txt'
	
	numberOfMeasurements = 0 #corresponding to the number of chromosome arms that we have, this is set using the probability file. 
	chromosomes = []
	chromosomeArms = []
	positions = [] #store the chromosome names and position
	variantPositions = []
	noiseLevels = [] #store the noise levels that we wish to explore for the measurements. Currently this is just one level (Estimated from the reference sample)
	possibleArms = [] #store the chromosome arm names for convenience
	simulationProbabilityFile = 'simulationProbabilities.txt' #make all these things a separate settings file later
	numberOfSamples = 10
	
	permutationRounds = 1
	simulationRounds = 1
	
	usedVariantSlots = [] #Make sure that we do not use the same variant again in another subclone during evolution, store the used positions here. 
	simulationProbabilities = None
	
	def initialize(self):
		#Pre-processing:
		#- read the probabilities of certain events (convert this to a better file format)
		simulationProbabilities = SimulationProbabilities()
		#- store the probabilities in an object (simulationprobabilities)
		simulationProbabilities.readFullProbabilityFile(self.simulationProbabilityFile)
		
		for armProbabilityList in simulationProbabilities.armProbabilities:
			self.possibleArms.append(armProbabilityList[0])

		#Parse the reference file so that we only have to read this one once.
		#We store the positions and chromosomes
		self.parseReferenceFile()
		
		#For each of the somatic variants, we know their positions.
		#We also know the positions of the measurements
		#Read the measurement positions and keep an index of where the somatic variant would go.
		#These values are the somatic variant indices. 
		self.variantPositions = self.obtainSomaticVariantIndices(simulationProbabilities)
		self.numberOfMeasurements = len(self.positions) #this is the total number of arms that we have in our test set.
		self.simulationProbabilities = simulationProbabilities

	
	#Parse the reference file and store the chromosomes/positions. Annotate a chromosome with p/q arm!
	
	def obtainSomaticVariantIndices(self, simulationProbabilities):
				
		offset = 0
		variantIndices = []
		for variant in simulationProbabilities.somaticVariants:
			
			position = variant.position
			chromosome = variant.chromosome

			for measurementPosition in range(0, len(self.positions)-1):
				
				if str(chromosome) == str(self.chromosomeArms[measurementPosition]): #the chromosome needs to be the same
					
					#for all the variants within the SNP range		
					if int(position) > int(self.positions[measurementPosition]) and int(position) < int(self.positions[measurementPosition + 1]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
						offset += 1
					#Situation where the somatic variant comes before the SNP measurements
					if str(chromosome) != str(self.chromosomeArms[measurementPosition-1]) and int(position) <= int(self.positions[measurementPosition]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
						offset += 1
					#Situation where the somatic variant comes after the SNP measurements
					if str(chromosome) != str(self.chromosomeArms[measurementPosition+1]) and int(position) >= int(self.positions[measurementPosition]):
						variantIndex = measurementPosition + offset
						variantIndices.append(variantIndex)
						offset += 1
		
		return variantIndices
	
	def parseReferenceFile(self):
		#Read the segmentation, use this to annotate chromosome positions with p/q arm.
		
		segmentation = Segmentation()
		segmentation.setSegmentationFromFile('pq_segmentation.txt')
		
		#Read the reference file
		with open(self.referenceFile, "r") as inFile:
			referenceMeasurements = []
			lineInd = 0
			for line in inFile:
				#remove newlines and split file
				line = line.replace(",", ".") #make sure that commas are always dots
				#check the data in this line.
				#Split the line by tab
				line = line.strip('\r\n')
				splitLine = re.split("\t", line)
				
				chromosome = splitLine[0]
				position = splitLine[1]
				
				self.chromosomes.append(chromosome)
				self.positions.append(position)
				referenceMeasurements.append(float(splitLine[2]))
				chromosomeArmName = segmentation.getSegmentName(chromosome, position, position)
				self.chromosomeArms.append(chromosomeArmName)
				lineInd += 1
		
		#Store the measurement positions and chromosome names (maybe later also convert these to p/q arm?)
		npMeasurements = np.asarray(referenceMeasurements)
		self.noiseLevels.append(np.std(npMeasurements))
		
		#Compute the standard deviation for the measurements
	
	def subclonalExpansion(self):
		simulationProbabilities = self.simulationProbabilities
		#start simulating.
		
		
		#Add the chromosome information. Our measurements start and stop at specific chromosome positions.
		#Here we can just use the first position of a chromosome segment.
		
		chromosomes = []
		for arm in self.possibleArms:
			chromosomeName = re.findall("(\d+)", arm)
			
			if len(chromosomeName) < 1:
				chromosomeName = arm
			else:
				chromosomeName = chromosomeName[0]
			
			chromosomes.append(chromosomeName)
		
		cloneNumbers = []
		#Subclone has the properties: C, A (and allele identifiers) and somatic variants

		randomClone = random.sample(range(0,self.simulationRounds), 1)[0]
		
		allCScores = []
		allAScores = []
		allMuScores = []
		allTreeScores = []
		
		#Store the averages over 1000 permutations
		allPCScores = []
		allPAScores = []
		allPMuScores = []
		allPTreeScores = []
		
		cPValues = []
		aPValues = []
		muPValues = []
		treePValues = []
		
		#Pre-define the keys for the profiles to obtain a good sorting for our figures
		
		cProfile = self.defineCProfile()
		#aProfile = self.defineAProfile(EventDistances(self.kmin, self.kmax))
		aProfile = dict() #for these kmin and kmax, there are too many combinations to show. 
		muProfile = dict()
		#The mu profile is a bit much, there are too many numbers and the plot will be practically empty. 
		#muProfile = self.defineMuProfile(0,100)
		#print muProfile
		
		#Load the targetClone instance from a pickled file
				
		pkl_file = open('targetClone.pkl', 'rb')
		targetClone = pickle.load(pkl_file)
		pkl_file.close()
		
		now = datetime.datetime.now()
		newDir = 'Simulations/' + now.isoformat()
		os.mkdir(newDir)
		for iteration in range(0,self.simulationRounds):
			
			#Re-set all samples here each time, each sample is likely getting updated in the method (and I think the others are as well, this may be bad.)
			#1. Generate a healthy subclone
			healthySubclone = Subclone()
			healthySubclone.C = [C([2,3])]*(self.numberOfMeasurements)

			healthySubclone.name = 'Healthy'
			
			healthyAlleleObjects = []
			prevArm = self.chromosomeArms[0]
			alleleA = 'A' + str(uuid.uuid4())
			alleleB = 'B' + str(uuid.uuid4())
			print self.numberOfMeasurements
			for newAlleleInd in range(0, self.numberOfMeasurements):
				
				if self.chromosomeArms[newAlleleInd] != prevArm:
					alleleA = 'A' + str(uuid.uuid4())
					alleleB = 'B' + str(uuid.uuid4())
					
				newAlleles = Alleles(1,2)
				newAlleles.alleleIdentifiers = [alleleA, alleleB] #we add alleles A and B
				healthyAlleleObjects.append(newAlleles)
				prevArm = self.chromosomeArms[newAlleleInd]
			
			healthySubclone.A = healthyAlleleObjects
			
			healthySubclone.somaticVariants = simulationProbabilities.somaticVariants
			#Also initialize the somatic variants, we have a total of 36 which all start out with a value of 0
			
			#A healthy subclone is the same as a healthy cell, with the same C and A everywhere, which are 2 and AB.
			
			#2. Ask the healthy subclone to expand with duplication. We only make 1 child here.
			preGCNIS = healthySubclone.duplicateGenome(self)
			preGCNIS.name = 'Pre-GCNIS'
			
			
			#we first duplicate the genome and return this subclone.
			#this subclone should again expand.
			
			#3. Introduce the first couple of events creating the first malginant precursor
			GCNIS = preGCNIS.generateMalignancy(simulationProbabilities, self)
			GCNIS.name = 'GCNIS'

			somVarArray = []
			for variant in healthySubclone.somaticVariants:
				if variant.value == 1:
					somVarArray.append(1)
				else:
					somVarArray.append(0)
			
			#Make some sample objects for the current subclones
			healthySample = Sample(None, None)
			healthySample.C = healthySubclone.C
			healthySample.A = healthySubclone.A
			healthySample.Mu = Mu(0)
			#obtain the chromosome, start and end information from the other samples
			measurements = self.generateMeasurements(healthySubclone, 100) #the mu of the tumor is here 0 because we have a healthy sample
			healthySample.afMeasurements = measurements[0]
			healthySample.measurements = measurements[1]
	
			healthySample.somaticVariants = somVarArray
			healthySample.somaticVariantsInd = self.variantPositions
			healthySample.setParent(None)
			healthySample.name = healthySubclone.name
			
			#Do we still need these values? Add them anyway for now
			eventDistances = EventDistances(self.kmin, self.kmax)
			bestCMuHealthy = CMuCombination(C([2,3]), Mu(0), eventDistances)
			
			healthySample.bestCMu = [bestCMuHealthy]*len(measurements[1].measurements)
			healthySample.originalCMu = healthySample.bestCMu
			
			####Pre-GCNIS
			#Also repeat this for pre-gcnis and gcnis, but then generate the LAF and mu randomly
			preGCNISSample = Sample(None, None)
			preGCNISSample.C = preGCNIS.C
			preGCNISSample.A = preGCNIS.A
			
			#Define a 4N precursor
			precursor4N = Subclone()
			precursor4N.C = [C([2,3])]*(self.numberOfMeasurements)

			precursor4N.name = 'Healthy'
			
			healthyAlleleObjects = []
			prevArm = self.chromosomeArms[0]
			alleleA = 'A' + str(uuid.uuid4())
			alleleB = 'B' + str(uuid.uuid4())
			for newAlleleInd in range(0, self.numberOfMeasurements):
				
				if self.chromosomeArms[newAlleleInd] != prevArm:
					alleleA = 'A' + str(uuid.uuid4())
					alleleB = 'B' + str(uuid.uuid4())
					
				newAlleles = Alleles(1,2)
				newAlleles.alleleIdentifiers = [alleleA, alleleB] #we add alleles A and B
				healthyAlleleObjects.append(newAlleles)
				prevArm = self.chromosomeArms[newAlleleInd]
			
			precursor4N.A = healthyAlleleObjects
			
			precursor4N.somaticVariants = simulationProbabilities.somaticVariants
			precursor4NSample = Sample(None, None)
			precursor4NSample.C = precursor4N.C
			precursor4NSample.A = precursor4N.A
			precursor4NSample.Mu = Mu(0)
			#obtain the chromosome, start and end information from the other samples
			measurements = self.generateMeasurements(precursor4N, 0) #the mu of the tumor is here 0 because we have a healthy sample
			precursor4NSample.afMeasurements = measurements[0]
			precursor4NSample.measurements = measurements[1]
	
			precursor4NSample.somaticVariants = somVarArray
			precursor4NSample.somaticVariantsInd = self.variantPositions
			precursor4NSample.setParent(None)
			precursor4NSample.name = precursor4N.name
		
			bestCMu4N = CMuCombination(C([2,3]), Mu(0), eventDistances)
			
			precursor4NSample.bestCMu = [bestCMu4N]*len(measurements[1].measurements)
			precursor4NSample.originalCMu = precursor4NSample.bestCMu
			
			
			somVarArray = []
			for variant in preGCNIS.somaticVariants:
				if variant.value == 1:
					somVarArray.append(1)
				else:
					somVarArray.append(0)
			
			preGCNISSample.somaticVariants = somVarArray
			preGCNISSample.somaticVariantsInd = self.variantPositions
			#Calculate these measurements from the counts
			measurements = self.generateMeasurements(preGCNIS, 100) #assume that this is the last precursor we know, which is then 100% tumor. 
			preGCNISSample.afMeasurements = measurements[0]
			preGCNISSample.measurements = measurements[1]
	
			preGCNISSample.setParent(precursor4NSample) #use None so that it will use itself as parent (first known subclone)
			preGCNISSample.name = 'Pre-GCNIS'
			randomMuT = random.sample(range(0, 100), 1)[0]#mu of the tumor component
			preGCNISSample.Mu =  Mu(0)
			
			#####GCNIS, repeat the same as above but then with GCNIS (would be nicer in a loop with a sample list but idc)
			GCNISSample = Sample(None, None)
			GCNISSample.C = deepcopy(GCNIS.C)
			GCNISSample.A = deepcopy(GCNIS.A)
			
			
			somVarArray = []
			for variant in GCNIS.somaticVariants:
				if variant.value == 1:
					somVarArray.append(1)
				else:
					somVarArray.append(0)
					
	
			GCNISSample.somaticVariants = somVarArray
			GCNISSample.somaticVariantsInd = self.variantPositions
			randomMuT = random.sample(range(0, 100), 1)[0]#mu of the tumor component
			#Calculate these measurements from the counts
			measurements = self.generateMeasurements(GCNIS, randomMuT)
			GCNISSample.afMeasurements = measurements[0]
			GCNISSample.measurements = measurements[1]
	
			GCNISSample.setParent(healthySample) #everyone's parent is the healthy sample in the first iteration because we do not have the information that we need yet! 
			GCNISSample.name = 'GCNIS'
			GCNISSample.Mu = Mu(0)
			
			allClones = [preGCNIS, GCNIS] #here we store all clones that we have made
			clones = self.evolve(GCNIS, 1, 1, simulationProbabilities, [])
			
			#select a random set of 10 clones
			
			#allClones = allClones + clones
			
			savedMu = [precursor4NSample.Mu, preGCNISSample.Mu, GCNISSample.Mu]
			cloneInd = 0
			
			#We want to select 10 samples from all other samples except for the precursor.
			
			samples = [precursor4NSample] #reset the samples each iteration
			precursorSamples = [preGCNISSample, GCNISSample]
			#1. select if we want GCNIS and preGCNIS or not.
			
			finalClones = [precursor4N]
			#numberOfPrecursorsToSample = random.sample(range(1, len(samples + precursorSamples)), 1)[0]
			numberOfPrecursorsToSample = 2 #no random selection of samples, add everything
			if numberOfPrecursorsToSample == 0:
				samples = [precursor4NSample]
				finalClones = [precursor4N]
			#if we select 1, randomly choose between preGCNIS and GCNIS
			if numberOfPrecursorsToSample == 1:
				randomPrecursor = random.sample(range(0,1), 1)[0]
				samples.append(precursorSamples[randomPrecursor])
				finalClones.append(allClones[randomPrecursor])
			#if we select 2, select both
			if numberOfPrecursorsToSample == 2:
				samples += precursorSamples
				finalClones += allClones
			
			extraSamplesNeeded = self.numberOfSamples - len(samples)
			if len(clones) == 0:
				extraSamplesNeeded = 0
			elif extraSamplesNeeded > len(clones):
				extraSamplesNeeded = len(clones)-1 #never use more samples than clones 
			randomMuT = random.sample(range(0, 100), 1)[0]#mu of the tumor component
			#we should keep track of the original tree as well. Who is the parent of whom?
			#We have links to sample numbers, but the actual indices are not that clear. Can we compare objects?
			edgeList = [(0, precursor4NSample.name, preGCNISSample.name), (0, preGCNISSample.name, GCNISSample.name)]
			#sampleInd = random.sample(range(0, len(clones)), extraSamplesNeeded) #we select x more until we have a total of 10 samples.
			sampleInd = range(0,extraSamplesNeeded) #sample the first 10 samples, but keep a relationship between the samples. 
			cloneInd = 0
			sampleNames = []
			print "number of clones generated: ", len(clones)
			for clone in clones: #make a new sample for each clone
				#if cloneInd in sampleInd:
				finalClones.append(clone)
				newSample = Sample(None, None)
				#We cannot set C, A and mu because we do not know it yet
				print "new clone"
				somVarArray = []
				for variant in clone.somaticVariants:
					if variant.value == 1:
						somVarArray.append(1)
					else:
						somVarArray.append(0)
						
				newSample.somaticVariants = somVarArray
				newSample.somaticVariantsInd = self.variantPositions
				
				#update the A and C, 2 and 4 is not possible.
				possibilities = ['loss']
				prevChrom = ''
				
				for pos in range(0, len(clone.C)):
					
					dummyLaf = min(clone.A[pos].ACount, clone.A[pos].BCount) / float((clone.A[pos].ACount + clone.A[pos].BCount))
					
					if (clone.C[pos].c[1] == 2 or clone.C[pos].c[1] == 4) and dummyLaf == 0.5:
						
						if self.chromosomeArms[pos] != prevChrom: #for a new chromosome arm we choose a new event, use the same chromosome for all other arms. 
							#select a gain or loss
							choice = np.random.choice(possibilities, 1)
							
							if choice == 'loss':
								newC = clone.C[pos].c[1] - 1
								
								newA = Alleles(clone.A[pos].ACount - 1, clone.A[pos].BCount)
							if choice == 'gain':
								 newC = clone.C[pos].c[1] + 1
								 newA = Alleles(clone.A[pos].ACount, clone.A[pos].BCount + 1)
						clone.C[pos].c[1] = newC
						clone.A[pos] = newA
					
						prevChrom = self.chromosomeArms[pos]
					
				randomMuT = random.sample(range(0, 100), 1)[0] #mu of the tumor component
				
				newSampleMu = Mu(100-randomMuT)
				newSampleMu = Mu(0) #set the mu to always be no contamination. 
				savedMu.append(newSampleMu)
				#Calculate these measurements from the counts
				measurements = self.generateMeasurements(clone, randomMuT)
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
			
			
		
			#Count the number of unique vertices, then make a tree.
			vertices = []
			for edge in edgeList:
				parent = edge[1]
				child = edge[2]
				if parent not in vertices:
					vertices.append(parent)
				if child not in vertices:
					vertices.append(child)
			
			realTree = Graph(vertices, set(edgeList), edgeList)
			print "real edges: ", edgeList
			print "samples: ", vertices
			print "sample names: ", sampleNames
			
			#For each of the clones, we need to start making samples.
			#A sample has:  - a name (this is random apart from for the precursor, pre-gcnis and gcnis)
			#				- afMeasurements
			#				- somatic variants
			#				- pre-defined mu
			#We could write the data to files, but for multiple iterations this may take a bit long.
		
			
			#obtain the C for every clone and merge it into a matrix
			
			cMatrix = np.empty([self.numberOfMeasurements, len(finalClones)], dtype=float)
			for cloneInd in range(0, len(finalClones)):
				#make the C from the individual Cs
				cloneC = finalClones[cloneInd].C
				for cInd in range(0, len(cloneC)):
					cMatrix[cInd][cloneInd] = cloneC[cInd].c[1]
			
			#Also print A
					
			aMatrix = np.empty([self.numberOfMeasurements, len(finalClones)], dtype=object)
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
						
			measurementsMatrix = np.empty([self.numberOfMeasurements, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].afMeasurements
				measurementsMatrix[:,cloneInd] = measurements
				
			lafMeasurementsMatrix = np.empty([self.numberOfMeasurements, len(samples)], dtype=float)
			for cloneInd in range(0, len(samples)):
				measurements = samples[cloneInd].measurements.measurements
				lafMeasurementsMatrix[:,cloneInd] = measurements
			
			#Somatic variants matrix
			np.savetxt(newDir + '/SomVar_' + str(iteration) + '.txt', somVarMatrix, fmt='%i', delimiter='\t')
			
			
			
		
			[eCMatrix, eAMatrix, eSamples, trees] = targetClone.run(samples)
			
			#1. Implement different scoring for the output
			
			eCMatrixFloat = eCMatrix.astype(float)
			cMatrixFloat = cMatrix.astype(float)
			#cScore = self.computeCRMSE(cMatrixFloat, eCMatrixFloat)
			cScore = self.computeCError(cMatrixFloat, eCMatrixFloat)
			
			allCScores.append(cScore / cMatrixFloat.size)
			
			##Also go through the C matrices and obtain the differences
			cProfile = self.determineCProfile(cProfile, eCMatrixFloat, cMatrixFloat)
			
			#aScore = self.computeARMSE(aMatrix, eAMatrix)
			aData = self.computeAError(aMatrix, eAMatrix)
			aScore = aData[0] / float(aMatrix.size)
		
			allAScores.append(aScore)
			eAStringMatrix = aData[2]
			aStringMatrix = aData[1]

			#Also obtain the profile for the alleles
			aProfile = self.determineAProfile(aProfile, eAMatrix, aMatrix)
			
			muData = self.computeMuError(eSamples, savedMu)
			muScore = muData[0]
			allMuScores.append(muScore)
			allRealMuT = muData[1]
			allEMuT = muData[2]
			
			#Compute the mu profile
			muProfile = self.determineMuProfile(muProfile, muData[2], muData[1])
			
			treeScore = self.computeTreeError(trees, realTree)
			allTreeScores.append(treeScore)
			bestTree = trees[len(trees)-1]
			#for each iteration we can save the individual scores. We use these later to make some plots of the current performance.
			
			#Also save the simulated data to files, we can then look back at these later if we need to have more information.
			
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
		
			#we should actually concatenate these throughout the iterations
			# max = 6
			
			# confMat = self.makeCConfusionMatrix(cMatrixFloat, eCMatrixFloat, self.kmin, max) #the max is actually 6, this is related to amplicons. 
			# self.plotConfusionMatrix(confMat, 'Output/confMat_test.png', "".join([str(i) for i in range(self.kmin, max+1)])) #does the alphabet need to be strings?
			
			#Compute the amibguity score. How often do we expect to be able to not solve the problem? What is our maximum expected score? 
			
			#For every real A, compute the real LAF in that position (without noise)
			#Compute the expected LAF for the inferred A as well
			#If there is an ambiguity, count that we could not solve the problem.
			
			[ambiguityScore, totalScore] = self.computeAmbiguityScore(aMatrix, eAMatrix, allRealMuT, allEMuT)
			
			print "ambiguity score: ", ambiguityScore
			
			
			print "total score without ambiguities: ", totalScore
			
			
			print "all C errors: ", allCScores
			
			
			print "all A errors: ", allAScores
			
			
			print "all mu errors: ", allMuScores
			
			
			print "all Tree errors: ", allTreeScores
			
			
		
		
			##Now that we have the errors of the method outcome on the actual simulated dataset, let's make permutations and compute the p-value.
			
			#1. Obtain a matrix with all the AF measurements of the samples
			#2. Shuffle the AF measurements within and between samples
			#3. Convert the AF measurements to LAF measurements and store the new information in the permuted sample objects. 
			#4. Obtain a matrix with all the somatic variants of the samples
			#5. Shuffle the somatic variants in the same way as the AF data
			#6. Update the samples
			#7. Run TargetClone and compute the RMSE for the permutation
		# 	
		# 	start_time = time.time()
		# 	#1. AF measurements
		# 	afMeasurementsMatrix = np.empty([len(samples[0].afMeasurements), len(samples)], dtype=float)
		# 	somaticVariantMatrix = np.empty([len(samples[0].somaticVariants), len(samples)], dtype=float)
		# 	for sampleInd in range(0, len(samples)):
		# 		
		# 		afMeasurementsMatrix[:,sampleInd] = deepcopy(samples[sampleInd].afMeasurements)
		# 		somaticVariantMatrix[:,sampleInd] = deepcopy(samples[sampleInd].somaticVariants)
		# 	
		# 	start_timeP = time.time()
		# 	#result_list = pool.map(process_x, my_list) #make permutations in a function, call this function here with multiple arguments
		# 	 #turn off permutations for now. 
		# 	
		# 	#return the object containing the scores. We can access those here. 
		# 	permutationObj = self.permutations(afMeasurementsMatrix, somaticVariantMatrix, samples)
		# 	
		# 	pCScores = []
		# 	pAScores = []
		# 	pMuScores = []
		# 	pTreeScores = []
		# 	#Count how many permutation statistics are smaller than the statistic. 
		# 	smallerCountC = 0
		# 	smallerCountA = 0
		# 	smallerCountMu = 0
		# 	smallerCountTree = 0
		# 	for permutationRound in range(0, len(permutationObj.pC)):
		# 		
		# 		pCMatrix = permutationObj.pC[permutationRound]
		# 		pAMatrix = permutationObj.pA[permutationRound]
		# 		pSamples = permutationObj.pSamples[permutationRound]
		# 		pTrees = permutationObj.pTrees[permutationRound]
		# 		
		# 		# np.savetxt(newDir + '/PC_' + str(permutationRound) + '.txt', cMatrix.astype(int), fmt='%i', delimiter='\t')
		# 		# 
		# 		# #Save A
		# 		# np.savetxt(newDir + '/PA_' + str(permutationRound) + '.txt', aStringMatrix, fmt='%s', delimiter='\t')
		# 		# 
		# 		# #Save Mu
		# 		# np.savetxt(newDir + '/PMu_' + str(permutationRound) + '.txt', allRealMuT, fmt='%f', delimiter='\t')
		# 		# 
		# 		# #Save the trees
		# 		# f = open(newDir + '/PTrees_' + str(permutationRound) + '.txt', 'w')
		# 		# f.write(str(realTree.getGraph()))  # python will convert \n to os.linesep
		# 		# f.close()
		# 		# 
		# 		pCMatrixFloat = pCMatrix.astype(float)
		# 		pCError = self.computeCError(cMatrixFloat, pCMatrixFloat)
		# 		
		# 		pCScores.append(pCError / pCMatrixFloat.size)
		# 		
		# 		aData = self.computeAError(aMatrix, pAMatrix)
		# 		pAError = aData[0] / float(aMatrix.size)
		# 		pAScores.append(pAError)
		# 		
		# 		muData = self.computeMuError(pSamples, savedMu)
		# 		pMuError = muData[0]
		# 		pMuScores.append(pMuError)
		# 		
		# 		pTreeScore = self.computeTreeError(pTrees, realTree)
		# 		pTreeScores.append(pTreeScore)
		# 		
		# 		#Check if the statistic is smaller than the observed statistic.
		# 		if pCError < cScore:
		# 			smallerCountC += 1
		# 		if pAError < aScore:
		# 			smallerCountA += 1
		# 		if pMuError < muScore:
		# 			smallerCountMu += 1
		# 		if pTreeScore < treeScore:
		# 			smallerCountTree += 1
		# 	#Compute the permutation error scores
		# 
		# 	print "permuted C errors: ", sum(pCScores) / len(pCScores)
		# 	print "permuted A errors: ", sum(pAScores) / len(pAScores)
		# 	print "permuted mu errors: ", sum(pMuScores) / len(pMuScores)
		# 	print "permuted tree errors: ", sum(pTreeScores) / len(pTreeScores)
		# 	[ambiguityScore, totalScore] = self.computeAmbiguityScore(aMatrix, pAMatrix, muData[1], muData[2])
		# 	
		# 	print "ambiguity score in permutations: ", ambiguityScore
		# 	
		# 	
		# 	print "total score without ambiguities in permutations: ", totalScore
		# 	
		# 	allPCScores.append(pCScores)
		# 	allPAScores.append(pAScores)
		# 	allPMuScores.append(pMuScores)
		# 	allPTreeScores.append(pTreeScores)
		# 	
		# 	# allPCScores.append(sum(pCScores) / len(pCScores))
		# 	# allPAScores.append(sum(pAScores) / len(pAScores))
		# 	# allPMuScores.append(sum(pMuScores) / len(pMuScores))
		# 	# allPTreeScores.append(sum(pTreeScores) / len(pTreeScores))
		# 	
		# 	#Compute some p-values as well
		# 	# cPValues.append(smallerCountC / float(self.permutationRounds))
		# 	# aPValues.append(smallerCountA / float(self.permutationRounds))
		# 	# muPValues.append(smallerCountMu / float(self.permutationRounds))
		# 	# treePValues.append(smallerCountTree / float(self.permutationRounds))
		# #Write the permutation scores to files. 
		# f = open(newDir + '/pCError.txt', 'w')
		# f.write(','.join([str(i) for i in allPCScores]))  # python will convert \n to os.linesep
		# f.close()
		# f = open(newDir + '/pAError.txt', 'w')
		# f.write(','.join([str(i) for i in allPAScores]))  # python will convert \n to os.linesep
		# f.close()
		# f = open(newDir + '/pMuError.txt', 'w')
		# f.write(','.join([str(i) for i in allPMuScores]))  # python will convert \n to os.linesep
		# f.close()
		# f = open(newDir + '/pTreeError.txt', 'w')
		# f.write(','.join([str(i) for i in allPTreeScores]))  # python will convert \n to os.linesep
		#f.close()
		
		#Write the error to a file. In the end, we wish to obtain 100 files, 1 for all the simulation errors, and 1 for all the permutation errors.
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
		
		
			
		#Show the profiles
		# print cProfile
		# print aProfile
		# print muProfile
		# self.printProfile(cProfile, 'cSignature.png', 'Signature of errors in C', 'Errors in C', 'Frequency')
		# #Repeat for A and mu
		# self.printProfile(aProfile, 'aSignature.png', 'Signature of errors in A', 'Errors in A', 'Frequency')
		# self.printProfile(muProfile, 'muSignature.png', 'Signature of errors in mu', 'Errors in mu', 'Frequency')
		# 
		#Make some plots of the scores! Run the method for asround 20 iterations and see what happens.
		# print "all C errors: ", allCScores
		# cData = [allCScores, allPCScores]
		# plt.figure()
		# plt.boxplot(cData)
		# #Add a title, axis labels
		# ax = axes()
		# ax.set_xticklabels(['Simulated C', 'Permuted C'])
		# plt.xlabel('')
		# plt.ylabel('Error')
		# plt.title('Error of the C matrix in 100 simulation datasets with 1000 permutations')
		# plt.savefig('cScores.png')
		# plt.clf()
		# 
		# print "all A errors: ", allAScores
		# aData = [allAScores, allPAScores]
		# plt.figure()
		# plt.boxplot(aData)
		# #Add a title, axis labels
		# ax = axes()
		# ax.set_xticklabels(['Simulated A', 'Permuted A'])
		# plt.xlabel('')
		# plt.ylabel('Error')
		# plt.title('Error of the A matrix in 100 simulation datasets with 1000 permutations')
		# plt.savefig('aScores.png')
		# plt.clf()
		# 
		# print "all mu errors: ", allMuScores
		# muData = [allMuScores, allPMuScores]
		# plt.figure()
		# plt.boxplot(muData)
		# #Add a title, axis labels
		# ax = axes()
		# ax.set_xticklabels(['Simulated Mu', 'Permuted Mu'])
		# plt.xlabel('')
		# plt.ylabel('Error')
		# plt.title('Error of the Mu matrix in 100 simulation datasets with 1000 permutations')
		# plt.savefig('muScores.png')
		# plt.clf()
		# 
		# print "all Tree errors: ", allTreeScores
		# treeData = [allTreeScores, allPTreeScores]
		# plt.figure()
		# plt.boxplot(treeData)
		# #Add a title, axis labels
		# ax = axes()
		# ax.set_xticklabels(['Simulated Trees', 'Permuted Trees'])
		# plt.xlabel('')
		# plt.ylabel('Edit distance')
		# plt.title('Tree edit distances in 100 simulation datasets with 1000 permutations')
		# plt.savefig('treeScores.png')
		# plt.clf()
		# 
		# print "p-values for C: ", cPValues
		# print "p-values for A: ", aPValues
		# print "p-values for Mu: ", muPValues
		# print "p-values for Trees: ", treePValues
		# 
		#print time.time() - start_time, " seconds to run one simulation with ", self.permutationRounds, " permutations"
		
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
	
	def permutations(self, afMatrix, somVarMatrix, samples):
				
		#Test the permutations
		pool = Pool(processes=8)

		permutationHandler = Permute(afMatrix, somVarMatrix, samples, self) #pass the simulation functions
		for permutation in range(0,self.permutationRounds):
			pool.apply_async(permutationHandler.permute())
			
		#print the scores to test the outcome
		return permutationHandler
	
	def makeCConfusionMatrix(self, cMatrixFloat, eCMatrixFloat, min, max):
		
		#Obtain all the unique values in both matrices
		kRange = range(min, max+1)
		confMat = np.zeros([len(kRange), len(kRange)])
		
		#Loop through the positions. Check which errors are commonly made.
		for row in range(0, cMatrixFloat.shape[0]):
			for col in range(0, cMatrixFloat.shape[1]):
				realC = int(cMatrixFloat[row][col])
				estimatedC = int(eCMatrixFloat[row][col])
				
				#Find the index of the values in the kRange
				#add it to the confusion matrix.
				realCInd = kRange.index(realC)
				estimatedCInd = kRange.index(estimatedC)
				#confMat[realCInd][estimatedCInd] += 1
				confMat[estimatedCInd][realCInd] += 1
		return confMat
	
	def plotConfusionMatrix(self, conf_arr, output, alphabet):
		norm_conf = []
		for i in reversed(conf_arr): #reverse it to switch around the axes of the plot
			a = 0
			tmp_arr = []
			a = sum(i, 0)
			for j in i:
				tmp_arr.append(j)
			norm_conf.append(tmp_arr)

		fig = plt.figure()
		plt.clf()
		ax = fig.add_subplot(111)
		ax.set_aspect(1)
		res = ax.imshow(np.array(norm_conf), cmap=plt.cm.jet, 
						interpolation='nearest')
		
		width, height = conf_arr.shape
		
		for x in xrange(width):
			for y in xrange(height):
				ax.annotate(str(norm_conf[x][y]), xy=(y, x), 
							horizontalalignment='center',
							verticalalignment='center')
		
		cb = fig.colorbar(res)
		plt.xticks(range(width), alphabet[:width])
		plt.yticks(range(height), list(reversed(alphabet))[:height])
		plt.savefig(output, format='png')
	
	def computeCError(self, cMatrixFloat, eCMatrixFloat):

		#Compute the absolute distance between the matrices.
		absDifference = np.absolute(eCMatrixFloat - cMatrixFloat)
		#Check which values are larger than 0.
		incorrectPositions = np.zeros(absDifference.shape)
		incorrectPositions[np.where(absDifference > 0)] = 1
		totalError = np.sum(incorrectPositions)
		return totalError
	
	def computeAError(self, aMatrix, eAMatrix):
		
		#Compute the edit distance between alleles. If the alleles are not the same, report an error.
		aStringMatrix = np.empty(aMatrix.shape, dtype=object)
		eAStringMatrix = np.empty(eAMatrix.shape, dtype=object)
		
		totalDifference = 0
		for cloneInd in range(0, aMatrix.shape[1]):
			for alleleInd in range(0,aMatrix.shape[0]):
				
				if aMatrix[alleleInd][cloneInd].ACount > aMatrix[alleleInd][cloneInd].BCount:
					oldACount = aMatrix[alleleInd][cloneInd].ACount
					oldBCount = aMatrix[alleleInd][cloneInd].BCount
					aMatrix[alleleInd][cloneInd].BCount = oldACount
					aMatrix[alleleInd][cloneInd].ACount = oldBCount
					
				aCount = aMatrix[alleleInd][cloneInd].ACount
				bCount = aMatrix[alleleInd][cloneInd].BCount

				aStringMatrix[alleleInd][cloneInd] = aMatrix[alleleInd][cloneInd].getAllelesAsString()
				if eAMatrix[alleleInd][cloneInd] is np.nan:
					continue
				
				eACount = eAMatrix[alleleInd][cloneInd].ACount
				eBCount = eAMatrix[alleleInd][cloneInd].BCount
				eAStringMatrix[alleleInd][cloneInd] = eAMatrix[alleleInd][cloneInd].getAllelesAsString()
				
				#Switch the alleles around if we have 'A' instead of 'B'.
				#So 'AAB' should be 'ABB'
				
				
				if abs(aCount - eACount) > 0 or abs(bCount - eBCount) > 0:
					totalDifference += 1
				
				#difference = abs(aCount - eACount) + abs(bCount - eBCount)
				#totalDifference += difference
		#aScore = np.sqrt(np.power((totalDifference / float(aMatrix.shape[0])), 2))
		
		return [totalDifference, aStringMatrix, eAStringMatrix]
	
	def computeMuError(self, eSamples, savedMu):
		allRealMuT = []
		allEMuT = []
		
		#Also check the difference between mu
		totalMuDifference = 0
		for sample in range(0, len(eSamples)):
			estimatedMuT = eSamples[sample].Mu.mu[1]
			realMuT = savedMu[sample].mu[1]
			
			muDifference = np.absolute(realMuT - estimatedMuT)
			
			allRealMuT.append(realMuT)
			allEMuT.append(estimatedMuT)
			#if muDifference > 0.5:
			totalMuDifference += muDifference
		
		muError = totalMuDifference / len(eSamples)
		
		return [muError, allRealMuT, allEMuT]
	
	def computeTreeError(self, trees, realTree):
		bestTree = trees[len(trees)-1]
		print "the best tree: ", bestTree.edgeList
		#compare the best tree to the edges of our tree.
		treeDistance = bestTree.computeTreeEditDistance(realTree)
		
		return treeDistance #use this for now. 
		

	def computeCRMSE(self, cMatrixFloat, eCMatrixFloat):
		#The current score shows how different the matrices are. 
		
		#Mask the NA positions
		#print eCMatrixFloat
		
		maskedECMatrixFloat = np.ma.array(eCMatrixFloat, mask = False)
		maskedECMatrixFloat.mask[np.isnan(eCMatrixFloat)] = True
		#print maskedECMatrixFloat
		
		cScore = np.sqrt(np.power((np.sum(np.absolute(maskedECMatrixFloat - cMatrixFloat)) / cMatrixFloat.size), 2))#we have to divide it by the sum of all values?
		return cScore
	
	def computeARMSE(self, aMatrix, eAMatrix):
		aStringMatrix = np.empty(aMatrix.shape, dtype=object)
		eAStringMatrix = np.empty(eAMatrix.shape, dtype=object)
		
		totalDifference = 0
		for cloneInd in range(0, aMatrix.shape[1]):
			for alleleInd in range(0,aMatrix.shape[0]):
				aCount = aMatrix[alleleInd][cloneInd].ACount
				bCount = aMatrix[alleleInd][cloneInd].BCount
				aStringMatrix[alleleInd][cloneInd] = aMatrix[alleleInd][cloneInd].alleleString
				if eAMatrix[alleleInd][cloneInd] is np.nan:
					continue
				
				eACount = eAMatrix[alleleInd][cloneInd].ACount
				eBCount = eAMatrix[alleleInd][cloneInd].BCount
				eAStringMatrix[alleleInd][cloneInd] = eAMatrix[alleleInd][cloneInd].alleleString
				difference = abs(aCount - eACount) + abs(bCount - eBCount)
				totalDifference += difference
		aScore = np.sqrt(np.power((totalDifference / float(aMatrix.shape[0])), 2))
		print aScore
		print "a score: ", aScore
		return aScore
	
	def computeMuRMSE(self, eSamples, savedMu):
		#Concatenate the mu as well so that we can obtain the profile for all mu at once
		allRealMuT = []
		allEMuT = []
		
		#Also check the difference between mu
		totalMuDifference = 0
		for sample in range(0, len(eSamples)):
			estimatedMuT = eSamples[sample].Mu.mu[1]
			realMuT = savedMu[sample].mu[1]
			
			print "estimated mu: ", estimatedMuT
			print "real mu T: ", realMuT
			
			totalMuDifference += np.absolute(realMuT - estimatedMuT)
			allRealMuT.append(realMuT)
			allEMuT.append(estimatedMuT)
		
		muScore = np.sqrt(np.power((totalMuDifference / 1), 2))
		print "mu score: ", muScore
		
		return [muScore, allRealMuT, allEMuT]
	
	def computeTreeRMSE(self, trees, realTree):
		#Finally, compare the trees. For now, we compare the distance to the final tree.
		bestTree = trees[len(trees)-1]
		print bestTree.getGraph()
		print realTree.getGraph()
		#compare the best tree to the edges of our tree.
		treeDistance = bestTree.computeTreeEditDistance(realTree)
		
		#Compute the score
		print "tree distance: ", treeDistance
		
		#Compute a RMSE score
		treeScore = np.sqrt(np.power((treeDistance / float(len(vertices))), 2))
		print "tree score: ", treeScore
		return treeScore
		
	def printProfile(self, profile, output, title, xlab, ylab):
		
		fig, ax = plt.subplots()
		fig.set_size_inches(18.5, 10.5)
		barWidth = 0.25
		width = 0.25
		opacity = 0.4

		cmap = plt.get_cmap('jet')
		colors = cmap(np.linspace(0, 1, 500)) #we don't know exactly how many colors we can use, set a higher number (12p issues?)
		
		previousVal = profile.keys()[0]
		previousVal = previousVal.split(",")[0]
		currentBar = []
		keyVals = []
		ticks = []
		offset = 0
		#For C, we wish to group the data by the k of the real C
		i = 1
		for profileEntry in profile.keys():
			#Make a new group each time the first value is within one group
			profileValues = profileEntry.split(",")
			firstVal = profileValues[0]
			if firstVal != previousVal: #append the frequency value to the current bar

				indices = [b + width + offset for b in np.arange(len(currentBar))]
				plt.bar(indices, currentBar, width, alpha=opacity, label=firstVal, color=colors[i]) #this may not be the nicest color to use
				#Save the ticks where we set specific values
				[ticks.append(index) for index in indices]
				offset += len(currentBar)
				#clear the values from the previous bar group
				currentBar = []
				keyVals = []
				barWidth += 0.25

				i += 1
			previousVal = firstVal
			currentBar.append(profile[profileEntry])
			keyVals.append(profileEntry)
		
		#Make sure that we also add the last one
		indices = [b + width + offset for b in np.arange(len(currentBar))]
		plt.bar(indices, currentBar, width, alpha=opacity, label=firstVal)
		#Save the ticks where we set specific values
		[ticks.append(index) for index in indices]
			
		#Also add the keys as the labels
		#the ticks are the total number of groups we made (i) * the number of bars.
		ax.set_xticks(ticks)
		ax.set_xticklabels(profile.keys(), rotation=45)
		ax.set_ylabel(ylab)
		ax.set_xlabel(xlab)
		ax.set_title(title)
		plt.tight_layout()
		plt.savefig(output, dpi=100)
		plt.clf()	
		
	
	def defineCProfile(self):
		#Go through all kmin and kmax
		#Make possible key combinations
		cProfile = OrderedDict()
		for k in range(self.kmin, self.kmax+1):
			for k2 in range(self.kmin, self.kmax+1):
				if k == k2:
					continue
				profileStr = str(k) + ',' + str(k2)
				cProfile[profileStr] = 0
		return cProfile
	
	def defineAProfile(self, eventDistances):
		#Make all possible combinations of alleles between kmin and kmax
		possibleAlleles = eventDistances.alleleList
		
		aProfile = OrderedDict()
		for a in possibleAlleles:
			for a2 in possibleAlleles:
				if a == a2:
					continue
				profileStr = a + ',' + a2
				aProfile[profileStr] = 0
		return aProfile
	
	def defineMuProfile(self, muMin, muMax):
		
		muProfile = OrderedDict()
		for mu in range(muMin, muMax):
			for mu2 in range(muMin, muMax):
				if mu == mu2:
					continue
				profileStr = str(mu/float(100)) + ',' + str(mu2/float(100))
				muProfile[profileStr] = 0
				
		return muProfile
	
	def determineCProfile(self, cProfile, eC, realC):
		#for each entry, check what the differences are
		for col in range(0, realC.shape[1]):
			for row in range(0, realC.shape[0]):
				
				realCVal = int(realC[row][col])
				if eC[row][col] is None or eC[row][col] is np.nan or math.isnan(eC[row][col]):
					eCVal = 'NaN'
				else:
					eCVal = int(eC[row][col])
					
				if realCVal == eCVal:
					continue
				profileStr = str(realCVal) + ',' + str(eCVal)
				
				#check if the value is already in the dictionary
				if profileStr not in cProfile.keys():
					cProfile[profileStr] = 0
				cProfile[profileStr] += 1 #incresae the frequency of this event. 
		return cProfile
	
	def determineAProfile(self, aProfile, eA, realA):
		
		#for each entry, check what the differences are
		for col in range(0, realA.shape[1]):
			for row in range(0, realA.shape[0]):
				
				#Obtain the allele string
				realAVal = realA[row][col].alleleString
				if eA[row][col] is None or eA[row][col] is np.nan:
					eAVal = 'NaN'
				else:
					eAVal = eA[row][col].alleleString
					
				if eAVal == realAVal:
					continue
					
				profileStr = str(realAVal) + ',' + str(eAVal)
				
				#check if the value is already in the dictionary
				if profileStr not in aProfile.keys():
					aProfile[profileStr] = 0
				aProfile[profileStr] += 1 #incresae the frequency of this event. 
		return aProfile
	
	def determineMuProfile(self, muProfile, eMuTList, realMuTList):
		
		#for each entry, check what the differences are
		for mu in range(0, len(realMuTList)):

			#Obtain the tumor mu values
			realMuT = realMuTList[mu]
			eMuT = eMuTList[mu]
			profileStr = str(realMuT) + ',' + str(eMuT)
			
			if realMuT == eMuT:
				continue
			
			#check if the value is already in the dictionary
			if profileStr not in muProfile.keys():
				muProfile[profileStr] = 0
			muProfile[profileStr] += 1 #incresae the frequency of this event. 
		return muProfile	
		
	####Here we end up duplicating our measurements, fix this!	
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
		#The LAF is easily computed for every position from the A matrix
		randomMuT = randomMuT / float(100)
		#The values here need to be weighed with a random mu, let's generate a random one here and see how it works out.
		#we also need to add the correct number of alleles for healthy cells to generate the correct AFs
		
		#Generate the measurement values from a normal distribution with added noise (first do this without noise)
		
		AF = []
		LAF = []
		#currentNoiseLevel = self.noiseLevels[0]
		
		currentNoiseLevel = 0
		#currentNoiseLevel = 0.02
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

			currentArm = self.chromosomeArms[measurementPosition] #at this point sample does not know where the measurements are as this is part of the measurement object
			armIndices = [i for i, x in enumerate(self.chromosomeArms) if x == currentArm]
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
			
			#Add noise to the measurements
			newAf = totalBCount / float(countSum) #we assume that the b count is the one that is the variant allele. We do not know!
			
			#take a sample from a normal distribution where the mu is the AF position (or should we sample from the LAF?)
			#noisedAf = np.random.normal(newAf, currentNoiseLevel, 1)[0]
			lower, upper = 0, 1 #AF truncated ones
			mu, sigma = newAf, currentNoiseLevel
			#X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
			#noisedAf = X.rvs(1)[0] #take a sample
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
		
		
	
	#Generate simulation data:
	#- Starting from a healthy cell
	#- Start from a healthy cell.
	#- The healthy cell duplicates into a 4N precursor
	#- The 4N precursor then undergoes binary subclonal expansion
	#- At each expansion, the subclones have a probability to aquire certain mutation.
	#- The possible mutations and their probabilities have been defined beforehand.
	#- If a mutation is not viable, the subclone 'dies' and is removed from further expansion.
	#- At the end we need to check how many subclones are left.
	#- When we have our subclones in place, we can start making sample objects of these and simulating measurements. 
	
	
class Subclone:
	
	C = None
	A = None
	somaticVariants = None
	parent = None
	children = None #we make a binary expansion, so there will always be only 1 new child, the other branch is the clone itself. 
	name = None
	
	def __init__(self):
		self.children = []
	
	def expand(self, duplicate):
		1+1
		
	def duplicateGenome(self, simulator):
		newC = []
		for c in self.C:
			oldCTumor = c.c[1]
			newC.append(C([2, oldCTumor*2]))
		newAlleles = []
		#Make sure that the allele identifiers are the same across a region with the same chromosome arm.
		
		prevArm = simulator.chromosomeArms[0]
		currentAlleleA = 'A' + str(uuid.uuid4())
		currentAlleleB = 'B' + str(uuid.uuid4())
		for a in range(0, len(self.A)):
			newACount = self.A[a].ACount*2
			newBCount = self.A[a].BCount*2
			#make the new names
			
			#Check for each arm if it is the same as the previous. If true, append the same name
			#if the chromosome arm is different, append a new one
			if simulator.chromosomeArms[a] != prevArm:
				currentAlleleA = 'A' + str(uuid.uuid4())
				currentAlleleB = 'B' + str(uuid.uuid4())
			
			#the identifiers for the new alleles depends on the ACount.
			newAlleleObjects = Alleles(newACount, newBCount)
			newAlleleObjects.alleleIdentifiers += self.A[a].alleleIdentifiers
			for newAIdentifier in range(0, (self.A[a].ACount*2 - self.A[a].ACount)):
				newAlleleObjects.alleleIdentifiers.append(currentAlleleA)
			for newBIdentifier in range(0, (self.A[a].BCount*2 - self.A[a].BCount)):
				newAlleleObjects.alleleIdentifiers.append(currentAlleleB)	
			newAlleles.append(newAlleleObjects)
			prevArm = simulator.chromosomeArms[a]
		
		#Check if our alleles always match across the chromosomes
		
		#we are using this for duplicating a precursor that does not have somatic variants. If the precursor did have somatic variants, we should duplicate these too here!
		
		precursor = Subclone()
		precursor.C = newC
		precursor.A = newAlleles
		precursor.somaticVariants = deepcopy(self.somaticVariants)
		precursor.parent = self
		precursor.name = str(uuid.uuid4())
		self.children.append(precursor)
		
		return precursor
		
	def generateMalignancy(self, simulationProbabilities, simulator):
		
		malignantSubclone = Subclone()
		malignantSubclone.C = deepcopy(self.C)
		malignantSubclone.A = deepcopy(self.A)
		malignantSubclone.somaticVariants = deepcopy(self.somaticVariants)
		malignantSubclone.parent = self
		self.children.append(malignantSubclone)
		malignantSubclone.name = str(uuid.uuid4())
		
		#Introduce 10 whole chromosome losses
		#From the probabilities we wish to obtain the possible chromosome numbers.
		#Draw 10 unique events. 
		#possibleArms = simulationProbabilities.armProbabilities.keys()
		possibleArms = []
		armLossProbabilities = []
		for armProbabilityList in simulationProbabilities.armProbabilities:
			possibleArms.append(armProbabilityList[0])
			armLossProbabilities.append(float(armProbabilityList[2]))
		
		#convert probabilities to numpy array	
		armLossProbabilities = np.asarray(armLossProbabilities)
		
		#These probabilities are not normalized, do this too:
		normalizedArmLossProbabilities = (armLossProbabilities / sum(armLossProbabilities))
		#print normalizedArmLossProbabilities
		
		chromosomes = []
		for arm in possibleArms:
			chromosomeName = re.findall("(\d+)", arm)
			
			if chromosomeName[0] == '12': #we skip 12
				continue 
			else:
				chromosomeName = chromosomeName[0]
			if chromosomeName not in chromosomes: #make sure that we do not add 1 for each arm but 1 for each chromosome8x
				chromosomes.append(chromosomeName)
		
		wholeChromosomesLost = random.sample(chromosomes, simulator.malignantWcLosses)

		#update the C and allele status of these chromosomes.
		#Each C entry corresponds to a chromosome arm. We need to link the chromosomes to the arms and lose those positions (1 copy at a time)

		previousChromosome = 0
		for chromosome in wholeChromosomesLost:
			for arm in range(0, len(possibleArms)):
				
				armChromosome = re.findall("(\d+)", possibleArms[arm])[0]
				
				
				if chromosome == str(armChromosome): #if this arm matches on the chromosome number

					chrInd = simulator.chromosomeArms.index(possibleArms[arm])
					chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[arm]]
					#Select the allele to be lost for the entire chromosome. All positions should behave the same!
					randomAllele = random.sample(malignantSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0] #select a random allele to be lost
					#select the allele from the list of identifiers. Here no somatic variants are present yet, so we can lose any copy. 
					randomAlelleIndex = malignantSubclone.A[chrIndices[0]].alleleIdentifiers.index(randomAllele)
	
					for ind in chrIndices:
						malignantSubclone.C[ind].c[1] = malignantSubclone.C[ind].c[1] -1

						del malignantSubclone.A[ind].alleleIdentifiers[randomAlelleIndex]
						#print "random allele: ", randomAllele
						#print "match 1: ", re.match('^A.+', randomAllele)
						#print "match 2: ", re.match('^B.+', randomAllele)
						#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount. 
						if re.match('^A.+', randomAllele) is not None:
							#print "pre-del A:", malignantSubclone.A[arm].ACount
							
							malignantSubclone.A[ind].ACount = malignantSubclone.A[ind].ACount - 1 #we always keep allele B
							#print "pro-del A:", malignantSubclone.A[arm].ACount
							#also remove the selected allele from the list
						if re.match('^B.+', randomAllele) is not None:
							
							malignantSubclone.A[ind].BCount = malignantSubclone.A[ind].BCount - 1 #we always keep allele A
							#print "pro-del B:", malignantSubclone.A[arm].ACount
					
		#10 arm losses
		#We select which arms are lost based on their probability: there is no probability compared to gain/loss, only relative to the other positions!
		
		#Here we should not select an arm that is completely lost (only here with malignancy, otherwise the clone in the expansion is not viable).
		#Check the current C counts, if there is one with a copy number of 1 then we can no longer lose it
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []

		for arm in range(0, len(possibleArms)):
			chrInd = simulator.chromosomeArms.index(possibleArms[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[arm]]
			
			if malignantSubclone.C[chrIndices[0]].c[1] > (simulator.kmin): #if there is only 1 copy left, we cannot lose any further. We select for this in the malignant precursor because it always needs to be viable!
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
		
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		armsLost = np.random.choice(filteredPossibleArms, simulator.cycleArmLosses, p=filteredNormalizedArmLossProbabilities, replace=False)
		#Update the lost arms in the C and A matrix
		for arm in range(0, len(armsLost)):
			
			#we need to find the right index, arms lost contains the values

			lostArmInd = possibleArms.index(armsLost[arm])
			
			chrInd = simulator.chromosomeArms.index(armsLost[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == armsLost[arm]]
			randomAllele = random.sample(malignantSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0] #select a random allele to be lost
			#select the allele from the list of identifiers. No somatic variants are introduced here yet, so we select any allele to be lost. 
			randomAlelleIndex = malignantSubclone.A[chrIndices[0]].alleleIdentifiers.index(randomAllele)

			for ind in chrIndices:
			
				malignantSubclone.C[ind].c[1] = malignantSubclone.C[ind].c[1] - 1

				#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount. 
				if re.match('^A.+', randomAllele) is not None:
					malignantSubclone.A[ind].ACount = deepcopy(malignantSubclone.A[ind].ACount) - 1 #we always keep allele B
					#also remove the selected allele from the list
				if re.match('^B.+', randomAllele) is not None:
					malignantSubclone.A[ind].BCount = deepcopy(malignantSubclone.A[ind].BCount) - 1 #we always keep allele A
				del malignantSubclone.A[ind].alleleIdentifiers[randomAlelleIndex]

					
		#Fill slots for somatic variants (starting from the previous precursor, then just update the correct objects)
		#20 random SVs
		
		somaticVariantsGained = random.sample(range(0, len(malignantSubclone.somaticVariants)), simulator.malignantSVGains)
		gainedPositions = []
		for var in somaticVariantsGained:
			gainedPositions.append(malignantSubclone.somaticVariants[var].ind)
		print "gained variants: ", gainedPositions
		simulator.usedVariantSlots += gainedPositions
		#Choose a random allele for each somatic variant, this is after duplication of the genome so we have many copies to choose from!
		for gainedVariant in somaticVariantsGained:
			#obtain the chromosome arm that these are located on
			chromosomeArm = simulationProbabilities.somaticVariants[gainedVariant].chromosome
			#Then select a random allele from this chromosome arm
			armIndex = possibleArms.index(chromosomeArm)
			chrInd = simulator.chromosomeArms.index(possibleArms[armIndex])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[armIndex]]
			
			randomAllele = random.sample(malignantSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0] #select a random allele to be lost
			
			#The somatic variant can be on any of these alleles.
			variantObject = malignantSubclone.somaticVariants[gainedVariant]
			variantObject.alleles.append(randomAllele)
			variantObject.value = 1
		
		#6 copies of chromosome 12p (we still need to get these values from the probability sheet)
		#Make sure that we also add allele identifiers!
		ch12pInd = simulator.chromosomeArms.index('12p')
		ch12pIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == '12p']
		randomAlleleInd = random.sample(range(0, len(malignantSubclone.A[ch12pInd].alleleIdentifiers)), 1)[0] #select a random allele to be lost, for 12p this is completely random.

		randomAllele = malignantSubclone.A[ch12pInd].alleleIdentifiers[randomAlleleInd]
		
		#Check if the random allele starts with A or B.
		#For either, set the prefix.
		
		newAlleleNames = []
		newAllelePrefix = 'A'
		if re.match('^A.+', randomAllele) is not None:
			
			newAllelePrefix = 'A'
			dupAlleleCopies = (6 - malignantSubclone.A[ch12pIndices[0]].BCount)
			
		if re.match('^B.+', randomAllele) is not None:
			
			newAllelePrefix = 'B'
			dupAlleleCopies = (6 - malignantSubclone.A[ch12pIndices[0]].ACount)
		#Set the identifiers of the new alleles
		for dupAlleleCopy in range(0, dupAlleleCopies): #add the extra copies of the duplicated allele. Give this allele a new identifier, otherwise we get issues with somatic variants
			newAlleleNames.append(newAllelePrefix + str(uuid.uuid4()))
	
		for ind in ch12pIndices:
			malignantSubclone.C[ind].c[1] = 6
		
			#We over-represent one of the alleles of 12p.
			if re.match('^A.+', randomAllele) is not None:
				malignantSubclone.A[ind].ACount = 6 - malignantSubclone.A[ind].BCount
				
				
			if re.match('^B.+', randomAllele) is not None:
				malignantSubclone.A[ind].BCount = 6 - malignantSubclone.A[ind].ACount
				
			#we have 6 copies in total, starting from 2 copies of A and two copies of B. The allele that is selected gets two extra copies with this identifier.
			
			malignantSubclone.A[ind].alleleIdentifiers += newAlleleNames
		
						
		# for a in range(0, len(malignantSubclone.A)):
		# 	print "chromosome arm: ", simulator.chromosomeArms[a]
		# 	print "alleles: ", malignantSubclone.A[a].alleleIdentifiers
		# exit()		
		return malignantSubclone
	
	def expandSubclone(self, simulationProbabilities, simulator): #simulation probabilities can be accessed through simulator, this is not necessary to provide!
		
		#Given the current subclone, we do a binary division
		#one subclone is the same as the current, the other one has mutations
		#we introduce:
		#- 8 losses based on their probability
		#- 3 gains based on their probability
		#- 2 somatic variants (ISA)
		#Always make sure that the kmin is 0 and the kmax is 5 and that these are not violated. (unless this is an amplification)
		
		#Do a test to see if the C are different or if we have references. 
		
		newSubclone = Subclone()
		newSubclone.C = deepcopy(self.C)
		newSubclone.A = deepcopy(self.A)
		newSubclone.somaticVariants = deepcopy(self.somaticVariants)
		newSubclone.parent = deepcopy(self)
		self.children.append(newSubclone)
		newSubclone.name = str(uuid.uuid4())
		
		
		##This part can go to the probability handler, we duplicate it now
		possibleArms = []
		armLossProbabilities = []
		for armProbabilityList in simulationProbabilities.armProbabilities:
			possibleArms.append(armProbabilityList[0])
			armLossProbabilities.append(float(armProbabilityList[2]))
		
		#convert probabilities to numpy array	
		armLossProbabilities = np.asarray(armLossProbabilities)
		
		#These probabilities are not normalized, do this too:
		normalizedArmLossProbabilities = (armLossProbabilities / sum(armLossProbabilities))
		
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []
		for arm in range(0, len(possibleArms)):
			chrInd = simulator.chromosomeArms.index(possibleArms[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[arm]]
			if newSubclone.C[chrIndices[0]].c[1] > (simulator.kmin-1): #we allow the loss of an arm, but if the total number of copies drops below kmin, the subclone is not viable.  
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
	
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		##Determine how many somatic variants we have per allele (the alleles should be the same across all indices of the chromosome arms)
		allelesAndVariantCounts = dict()
		for arm in possibleArms:
			armIndex = possibleArms.index(arm)
			chrInd = simulator.chromosomeArms.index(arm)
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == arm]
			for variant in newSubclone.somaticVariants:
				#Get the somatic variants that are present on this arm
				identifiers = variant.alleles
				if variant.chromosome == arm:
					#get the alleles
					
					for identifier in identifiers:
						if identifier not in allelesAndVariantCounts.keys():
							allelesAndVariantCounts[identifier] = []
							
						allelesAndVariantCounts[identifier].append(variant)
				else:
					for identifier in newSubclone.A[chrIndices[0]].alleleIdentifiers:
						allelesAndVariantCounts[identifier] = []
		
		probableElements = np.where(filteredNormalizedArmLossProbabilities > 0)[0]
		#print "pe: ", probableElements
		#Find how many elements in the normalized array are larger than 0. 
		armsToSample = 0
		if len(probableElements) < simulator.cycleArmLosses:
			armsToSample = len(probableElements)
		else:
			armsToSample = simulator.cycleArmLosses
		#print "to sample: ", armsToSample
		#print "length of possible arms: ", len(filteredPossibleArms)
		#print "length of normalized probs: ", filteredNormalizedArmLossProbabilities
		armsLost = np.random.choice(filteredPossibleArms, armsToSample, p=filteredNormalizedArmLossProbabilities, replace=False)

		#Update the lost arms in the C and A matrix
		for arm in range(0, len(armsLost)):
			#we need to find the right index, arms lost contains the values
			lostArmInd = possibleArms.index(armsLost[arm])
			chrInd = simulator.chromosomeArms.index(possibleArms[lostArmInd])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[lostArmInd]]
			
			#From all alleles, check which one has the least somatic variants. If there is a tie, we choose a random one.
			leastSVAllele = ''
			currentBest = (simulator.numberOfSVs + 1) #we can never has less than the maximum number of SVs
			#print "we have identifiers: ", newSubclone.A[lostArmInd].alleleIdentifiers
			#if there is no identifier, we have no alleles left and the clone is not viable.
			if len(newSubclone.A[chrIndices[0]].alleleIdentifiers) < 1:
				return False
			
			for allele in newSubclone.A[chrIndices[0]].alleleIdentifiers:
				#check how many somatic variants the alleles have
				if len(allelesAndVariantCounts[allele]) < currentBest:
					leastSVAllele = allele
				#randomAllele = random.sample(newSubclone.A[lostArmInd].alleleIdentifiers, 1)[0] #select a random allele to be lost
			
			#here we wish to remove the allele with the least number of somatic variants
				
			randomAlelleIndex = newSubclone.A[chrIndices[0]].alleleIdentifiers.index(leastSVAllele)
		
			for ind in chrIndices:
				
				newSubclone.C[ind].c[1] = newSubclone.C[ind].c[1] - 1 #only lose the copy after we are sure that we can lose the allele as well
				
				del newSubclone.A[ind].alleleIdentifiers[randomAlelleIndex]
				#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount.
				# print "prev allele A count: ", newSubclone.A[lostArmInd].ACount
				# print "prev allele B count: ", newSubclone.A[lostArmInd].BCount
				# print "current alleles: ", newSubclone.A[lostArmInd].alleleIdentifiers
				# print "losing allele: ", leastSVAllele
				# 
				if re.match('^A.+', leastSVAllele) is not None:
					newSubclone.A[ind].ACount = newSubclone.A[ind].ACount - 1 #we always keep allele B
					#also remove the selected allele from the list
				if re.match('^B.+', leastSVAllele) is not None:
					newSubclone.A[ind].BCount = newSubclone.A[ind].BCount - 1 #we always keep allele A
			
			lostVariants = allelesAndVariantCounts[leastSVAllele]
			
			varCounter = 0
			for lostVariant in lostVariants:
				print "lost variant: ", lostVariant.ind
				lostVariant.lost = True
				lostVariant.value = 0
				varCounter += 1
				
		#Gain alleles: it is always possible to gain, but only to gain alleles that have not yet been lost
		
		#The alleles to gain also depends on which alleles can still be gained.
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []
		for arm in range(0, len(possibleArms)):
			chrInd = simulator.chromosomeArms.index(possibleArms[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == possibleArms[arm]]
			if newSubclone.C[chrIndices[0]].c[1] < (simulator.kmax+1): #we allow the loss of an arm, but if the total number of copies drops below kmin, the subclone is not viable.  
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
		
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		#If we gain the allele that contains a somatic variant, we need to duplicate these as well. 
		armsGained = np.random.choice(filteredPossibleArms, simulator.cycleArmGains, p=filteredNormalizedArmLossProbabilities, replace=False)
		
		for arm in range(0, len(armsGained)):
			#The selected allele needs to be the same for all chromosomes. We only gain the allele once and always the same type!
			
			gainedArm = possibleArms.index(armsGained[arm])
			chrInd = simulator.chromosomeArms.index(armsGained[arm])
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == armsGained[arm]]
			
			aAlleles = []
			bAlleles = []

			for allele in newSubclone.A[chrIndices[0]].alleleIdentifiers:
				#if the arm matches on A, we add it to the A alleles, otherwise the B alleles.
				
				if re.match('^A.+', allele) is not None:
					aAlleles.append(allele)
					#print "a alleles: ", aAlleles
				else:
					bAlleles.append(allele)
					#print "b alleles: ", bAlleles
			
			#Here we do not allow an allele to be gained when the allele has already been lost!
			newAllelePrefix = 'A'
			if newSubclone.A[chrIndices[0]].ACount < 1 and newSubclone.A[chrIndices[0]].BCount > 0:
				for ind in chrIndices:
					newSubclone.A[ind].BCount = newSubclone.A[ind].BCount + 1
				#also choose which one to gain out of the B alleles
				randomAllele = random.sample(bAlleles, 1)[0]
				#then when we have selected one, we duplicate it in the allele identifiers
				newAllelePrefix = 'B'
			elif newSubclone.A[chrIndices[0]].BCount < 1 and newSubclone.A[chrIndices[0]].ACount > 0:
				for ind in chrIndices:
					newSubclone.A[ind].ACount = newSubclone.A[ind].ACount + 1
				randomAllele = random.sample(aAlleles, 1)[0]
			elif newSubclone.A[chrIndices[0]].BCount > 0 and newSubclone.A[chrIndices[0]].ACount > 0: #this should generally always be the case if we make sure that the kmin is not exceeded
				#select any old allele
				#why is there sometimes no allele here? Either of the above should work. This only means that sometimes even
				#print "identifier: ", newSubclone.A[gainedArm].alleleIdentifiers
				randomAllele = random.sample(newSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0]
				#print randomAllele
				if re.match('^A.+', randomAllele) is not None:
					for ind in chrIndices:
						newSubclone.A[ind].ACount = newSubclone.A[ind].ACount + 1 #we always keep allele B
				if re.match('^B.+', randomAllele) is not None:
					for ind in chrIndices:
						newSubclone.A[ind].BCount = newSubclone.A[ind].BCount + 1 #we always keep allele A
					newAllelePrefix = 'B'
			else: #if we have no copies left, we cannot do anything and w should skip this gain.
				continue
			newAlleleName = newAllelePrefix + str(uuid.uuid4())
			for ind in chrIndices:
				newSubclone.A[ind].alleleIdentifiers.append(newAlleleName)
			#we also notify the somatic variant that it has been duplicated
			duplicatedVariants = allelesAndVariantCounts[randomAllele]
			for duplicatedVariant in duplicatedVariants:
				duplicatedVariant.alleles.append(newAlleleName)
				
			#Only update C when we get here
			for ind in chrIndices:
				newSubclone.C[ind].c[1] = newSubclone.C[ind].c[1] + 1
		
		#Gain somatic variants: only gain somatic variants that have not been lost yet! We need to give these a tag to see if they have been lost or not in a previous iteration.
		#Currently losses are not implemented yet, so we can go with everything that is not already there.
		
		availableSomVarSlots = []
		for somVar in newSubclone.somaticVariants:
			#also check if the allele is still available, otherwise we cannot gain an allele there. 
			chromosome = somVar.chromosome
			chromosomeInd = possibleArms.index(somVar.chromosome)
			chrInd = simulator.chromosomeArms.index(somVar.chromosome)
			chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == somVar.chromosome]
			#print "the content: ", newSubclone.A[chromosomeInd].alleleIdentifiers
			if len(newSubclone.A[chrIndices[0]].alleleIdentifiers) < 1: #if there are no alleles to add the variant to, we skip it. 
				continue
			
			if somVar.value < 1 and somVar.ind not in simulator.usedVariantSlots: #Make sure that the position has not been used already.
				availableSomVarSlots.append(somVar)
		
		
		somaticVariantsGained = []
		if len(availableSomVarSlots) == 1: #if we cannot collect 2 variants to gain but at least one, we also need to sample only one
			somaticVariantsGained = random.sample(availableSomVarSlots, 1)
		if len(availableSomVarSlots) > 1:
			somaticVariantsGained = random.sample(availableSomVarSlots, simulator.cycleSVGains)
		
		if len(somaticVariantsGained) > 0: #if we can gain anything at all
			#somaticVariantsGained = random.sample(availableSomVarSlots, simulator.cycleSVGains)
			
			for somaticVariantGained in somaticVariantsGained:
				somaticVariantGained.value = 1
				simulator.usedVariantSlots.append(somaticVariantGained.ind)
				#it knows on which chromosome it has been gained on, but we still need to assign an allele identifier
				#if the chromosome arm already has alleles, we can use these
				#we consider the entire arm to consist of one big allele
				variantArm = somaticVariantGained.chromosome
				variantArmInd = possibleArms.index(variantArm)
				chrIndices = [i for i, x in enumerate(simulator.chromosomeArms) if x == variantArm]
				#Check the alleles that it has and select a random one
				randomAllele = random.sample(newSubclone.A[chrIndices[0]].alleleIdentifiers, 1)[0]
				somaticVariantGained.alleles.append(randomAllele)
		
		gainedPositions = []
		for var in somaticVariantsGained:
			gainedPositions.append(var.ind)
		print "gained variants in subclone: ", gainedPositions
		
		
		#Also check for viability, report a message to the method when we cannot go any further from here because the subclone is not viable (a bit of a recursive fashion).
		#Our viability checks are:
		#- We may not exceed our kmin and kmax
		#- 12p amplifications are not allowed to be lost (is this possible now?)
		i = 0
		for c in newSubclone.C:
			#print possibleArms[i]
			if simulator.chromosomeArms[i] == '12p': #exclude 12p from the viability check
				if c.c[1] > 6:
					return False #make sure that 12p is never gained further than 6 copies. 
				i += 1
				continue
			if c.c[1] < simulator.kmin or c.c[1] > simulator.kmax:
				return False
			i += 1
		
		# for a in range(0, len(newSubclone.A)):
		# 	print "chromosome arm: ", simulator.chromosomeArms[a]
		# 	print "alleles: ", newSubclone.A[a].alleleIdentifiers
		# exit()	
		return newSubclone
		
	
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
		

#file = 'lossGainProbabilityFile.txt'
#outfile = 'normalizedLossGainProbabilityFile.txt'
#SimulationProbabilities().convertProbabilities(file, outfile)

#SimulationProbabilities().readFullProbabilityFile('simulationProbabilities.txt')
Simulator()
	