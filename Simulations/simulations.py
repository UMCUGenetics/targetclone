#Classes for generating simulation data

import numpy as np
import csv
import re
import random
import uuid
import matplotlib.pyplot as plt
import scipy.stats as stats
from collections import OrderedDict
import math
import datetime
import os

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
	kmax = 5
	cellCycles = 6
	
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
	
	def __init__(self):
		#Pre-processing:
		#- read the probabilities of certain events (convert this to a better file format)
		simulationProbabilities = SimulationProbabilities()
		#- store the probabilities in an object (simulationprobabilities)
		simulationProbabilities.readFullProbabilityFile(self.simulationProbabilityFile)
		
		print simulationProbabilities.armProbabilities
		
		
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
		self.subclonalExpansion(simulationProbabilities)
	
	
	#Parse the reference file and store the chromosomes/positions. Annotate a chromosome with p/q arm!
	
	def obtainSomaticVariantIndices(self, simulationProbabilities):
				
		offset = 0
		variantIndices = []
		print "variants: ", len(simulationProbabilities.somaticVariants)
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
	
	def subclonalExpansion(self, simulationProbabilities):
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
		
		print chromosomes
		
		#Test the method first on two samples
		

		#add the healthy samples and gcnis samples to the c matrices and a matrices as well
		
		
		#For each clone, we now generate a sample object like the way we did above
		
		cloneNumbers = []
		#Subclone has the properties: C, A (and allele identifiers) and somatic variants
		iterations = 100
		randomClone = random.sample(range(0,iterations), 1)[0]
		
		allCScores = []
		allAScores = []
		allMuScores = []
		allTreeScores = []
		
		#Pre-define the keys for the profiles to obtain a good sorting for our figures
		
		cProfile = self.defineCProfile()
		print cProfile
		#aProfile = self.defineAProfile(EventDistances(self.kmin, self.kmax))
		aProfile = dict() #for these kmin and kmax, there are too many combinations to show. 
		print aProfile
		muProfile = dict()
		#The mu profile is a bit much, there are too many numbers and the plot will be practically empty. 
		#muProfile = self.defineMuProfile(0,100)
		#print muProfile
		
		
		targetClone = TargetClone() #pre-initialize the method's C and mu combinations, required only once
		now = datetime.datetime.now()
		newDir = 'Simulations/' + now.isoformat()
		os.mkdir(newDir)
		for iteration in range(0,iterations):
			
			#Re-set all samples here each time, each sample is likely getting updated in the method (and I think the others are as well, this may be bad.)
			#1. Generate a healthy subclone
			healthySubclone = Subclone()
			healthySubclone.C = [C([2,2])]*(self.numberOfMeasurements-2)
			#Separately handle chromosomes X and Y, these have different chromosome numbers and allele counts
			healthySubclone.C.append(C([2,1]))
			healthySubclone.C.append(C([2,1]))
			healthySubclone.name = 'Healthy'
			
			healthyAlleleObjects = []
			for newAlleleInd in range(0, self.numberOfMeasurements-2):
				newAlleles = Alleles(1,1)
				newAlleles.alleleIdentifiers = ['A' + str(uuid.uuid4()), 'B' + str(uuid.uuid4())] #we add alleles A and B
				healthyAlleleObjects.append(newAlleles)
			
			healthySubclone.A = healthyAlleleObjects
			xAlleles = Alleles(0,1)
			xAlleles.alleleIdentifiers = ['B' + str(uuid.uuid4())]
			healthySubclone.A.append(xAlleles)
			#make the object twice to ensure that we do not reference the other one
			yAlleles = Alleles(0,1)
			yAlleles.alleleIdentifiers = ['B' + str(uuid.uuid4())]
			healthySubclone.A.append(yAlleles)
			healthySubclone.somaticVariants = simulationProbabilities.somaticVariants
			#Also initialize the somatic variants, we have a total of 36 which all start out with a value of 0
			
			#A healthy subclone is the same as a healthy cell, with the same C and A everywhere, which are 2 and AB.
			
			#2. Ask the healthy subclone to expand with duplication. We only make 1 child here.
			preGCNIS = healthySubclone.duplicateGenome()
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
			healthySample.Mu = Mu(100)
			#obtain the chromosome, start and end information from the other samples
			measurements = self.generateMeasurements(healthySubclone, 0) #the mu of the tumor is here 0 because we have a healthy sample
			healthySample.afMeasurements = measurements[0]
			healthySample.measurements = measurements[1]
	
			healthySample.somaticVariants = somVarArray
			healthySample.somaticVariantsInd = self.variantPositions
			healthySample.setParent(None)
			healthySample.name = healthySubclone.name
			
			#Do we still need these values? Add them anyway for now
			eventDistances = EventDistances(self.kmin, self.kmax)
			bestCMuHealthy = CMuCombination(C([2,2]), Mu(100), eventDistances)
			
			healthySample.bestCMu = [bestCMuHealthy]*len(measurements[1].measurements)
			healthySample.originalCMu = healthySample.bestCMu
			
			####Pre-GCNIS
			#Also repeat this for pre-gcnis and gcnis, but then generate the LAF and mu randomly
			preGCNISSample = Sample(None, None)
			preGCNISSample.C = preGCNIS.C
			preGCNISSample.A = preGCNIS.A
			
			
			somVarArray = []
			for variant in preGCNIS.somaticVariants:
				if variant.value == 1:
					somVarArray.append(1)
				else:
					somVarArray.append(0)
			
			preGCNISSample.somaticVariants = somVarArray
			preGCNISSample.somaticVariantsInd = [1]*len(somVarArray)
			randomMuT = random.sample(range(0, 100), 1)[0] #mu of the tumor component
			#Calculate these measurements from the counts
			measurements = self.generateMeasurements(preGCNIS, randomMuT)
			preGCNISSample.afMeasurements = measurements[0]
			preGCNISSample.measurements = measurements[1]
	
			preGCNISSample.setParent(healthySample)
			preGCNISSample.name = 'Pre-GCNIS'
			preGCNISSample.Mu = Mu(100-randomMuT)
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
			GCNISSample.Mu = Mu(100-randomMuT)
			
			
			print "new method iteration: ", iteration
			allClones = [healthySubclone, preGCNIS, GCNIS] #here we store all clones that we have made
			clones = self.evolve(GCNIS, 1, 1, simulationProbabilities, [])
			
			allClones = allClones + clones
			
			savedMu = [healthySample.Mu, preGCNISSample.Mu, GCNISSample.Mu]
			cloneInd = 0
			samples = [healthySample, preGCNISSample, GCNISSample] #reset the samples each iteration
			
			randomMuT = random.sample(range(0, 100), 1)[0]#mu of the tumor component
			print "selected random mu: ", randomMuT
			#we should keep track of the original tree as well. Who is the parent of whom?
			#We have links to sample numbers, but the actual indices are not that clear. Can we compare objects?
			edgeList = [(0, healthySample.name, preGCNISSample.name), (0, preGCNISSample.name, GCNISSample.name)]
			for clone in clones: #make a new sample for each clone

				newSample = Sample(None, None)
				#We cannot set C, A and mu because we do not know it yet
				
				somVarArray = []
				for variant in clone.somaticVariants:
					if variant.value == 1:
						somVarArray.append(1)
					else:
						somVarArray.append(0)
						
				newSample.somaticVariants = somVarArray
				newSample.somaticVariantsInd = self.variantPositions
				
				randomMuT = random.sample(range(0, 100), 1)[0] #mu of the tumor component
				print "selected random mu: ", randomMuT
				newSampleMu = Mu(100-randomMuT)
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
				cloneInd += 1
			
			print edgeList

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
			#exit()
			
			#For each of the clones, we need to start making samples.
			#A sample has:  - a name (this is random apart from for the precursor, pre-gcnis and gcnis)
			#				- afMeasurements
			#				- somatic variants
			#				- pre-defined mu
			#We could write the data to files, but for multiple iterations this may take a bit long.
		
			
			#Use these data later to check if the method predictions match the actual values!
			cloneNumbers.append(len(clones))
			
			#obtain the C for every clone and merge it into a matrix
			
			cMatrix = np.empty([self.numberOfMeasurements, len(allClones)], dtype=float)
			for cloneInd in range(0, len(allClones)):
				#make the C from the individual Cs
				cloneC = allClones[cloneInd].C
				for cInd in range(0, len(cloneC)):
					cMatrix[cInd][cloneInd] = cloneC[cInd].c[1]
			
			#Also print A
					
			aMatrix = np.empty([self.numberOfMeasurements, len(allClones)], dtype=object)
			for cloneInd in range(0, len(allClones)):
				cloneA = allClones[cloneInd].A
				for aInd in range(0, len(cloneA)):
			
					aMatrix[aInd][cloneInd] = cloneA[aInd]
			
			#Check the somatic variants
			somVarMatrix = np.empty([len(allClones[1].somaticVariants), len(allClones)], dtype=float)
			for cloneInd in range(0, len(allClones)):
				somVar = allClones[cloneInd].somaticVariants
				for variant in range(0, len(somVar)):
					if somVar[variant].value == 1:
						somVarMatrix[variant][cloneInd] = 1
					else:
						somVarMatrix[variant][cloneInd] = 0
			#print somVarMatrix
			
			
			print "number of samples: ", len(samples)
			
			[eCMatrix, eAMatrix, eSamples, trees] = targetClone.run(samples)
			
			print "number of samples: ", len(samples)
			print "number of clones: ", len(allClones)
			
			print "number of returned samples: ", len(eSamples)
			print "size of real C: ", cMatrix.shape
			print "size of inferred C: ", eCMatrix.shape
				
			#If the method has introduced a precursor, we need to make sure that the matrices are compared with the same sizes.
			#if eCMatrix.shape[1] > cMatrix.shape[1]:
			#	eCMatrix = eCMatrix[:,1:cMatrix.shape[1]] #remove the last column for now, we only have 1 column extra
			
			#Compare the values to our real values.
			
			# print eCMatrix
			# print cMatrix
			# print cMatrix.size
			# 
			#The current score shows how different the matrices are. 
			eCMatrixFloat = eCMatrix.astype(float)
			cMatrixFloat = cMatrix.astype(float)
			#Mask the NA positions
			#print eCMatrixFloat
			
			maskedECMatrixFloat = np.ma.array(eCMatrixFloat, mask = False)
			maskedECMatrixFloat.mask[np.isnan(eCMatrixFloat)] = True
			#print maskedECMatrixFloat
			
			cScore = np.sqrt(np.power((np.sum(np.absolute(maskedECMatrixFloat - cMatrixFloat)) / cMatrix.size), 2))#we have to divide it by the sum of all values?
			print cScore
			print "c score: ", cScore
			
			allCScores.append(cScore)
			
			##Also go through the C matrices and obtain the differences
			cProfile = self.determineCProfile(cProfile, eCMatrixFloat, cMatrixFloat)
			
			
			#For the aMatrix we will need to compare the results a little bit different
			#We can instead sum the difference between the A and B counts (as in event distance?).
			
			#print aMatrix.shape
			#print eAMatrix.shape
			#print len(allClones)
			aStringMatrix = np.empty(aMatrix.shape, dtype=object)
			eAStringMatrix = np.empty(eAMatrix.shape, dtype=object)
			
			totalDifference = 0
			for cloneInd in range(0, len(allClones)):
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
			
			allAScores.append(aScore)
			
			#Also obtain the profile for the alleles
			aProfile = self.determineAProfile(aProfile, eAMatrix, aMatrix)
			
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
			allMuScores.append(muScore)
			
			#Compute the mu profile
			muProfile = self.determineMuProfile(muProfile, allEMuT, allRealMuT)
			
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
			allTreeScores.append(treeScore)
			#for each iteration we can save the individual scores. We use these later to make some plots of the current performance.
			
			#Also save the simulated data to files, we can then look back at these later if we need to have more information.
			
			
			
			
			#Save C
			np.savetxt(newDir + '/RealC_' + str(iteration) + '.txt', cMatrix.astype(int), fmt='%i', delimiter='\t')
			np.savetxt(newDir + '/EstimatedC_' + str(iteration) + '.txt', eCMatrixFloat, fmt='%f', delimiter='\t')
			
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
			
			
			##Now that we have the errors of the method outcome on the actual simulated dataset, let's make permutations and compute the p-value.
			
			

		#Show the profiles
		print cProfile
		print aProfile
		print muProfile
		self.printProfile(cProfile, 'cSignature.png', 'Signature of errors in C', 'Errors in C', 'Frequency')
		#Repeat for A and mu
		self.printProfile(aProfile, 'aSignature.png', 'Signature of errors in A', 'Errors in A', 'Frequency')
		self.printProfile(muProfile, 'muSignature.png', 'Signature of errors in mu', 'Errors in mu', 'Frequency')
		
		#Make some plots of the scores! Run the method for asround 20 iterations and see what happens.
		print "all C scores: ", allCScores
		plt.boxplot(allCScores)
		plt.xlabel('')
		plt.ylabel('RMSE')
		plt.title('RMSE of the C matrix in 100 simulation datasets')
		plt.savefig('cScores.png')
		plt.clf()
		print "all A scores: ", allAScores
		plt.boxplot(allAScores)
		plt.xlabel('')
		plt.ylabel('RMSE')
		plt.title('RMSE of the A matrix in 100 simulation datasets')
		plt.savefig('aScores.png')
		plt.clf()
		print "all mu scores: ", allMuScores
		plt.boxplot(allMuScores)
		plt.xlabel('')
		plt.ylabel('RMSE')
		plt.title('RMSE of the mu in 100 simulation datasets')
		plt.savefig('muScores.png')
		plt.clf()
		print "all Tree scores: ", allTreeScores
		plt.boxplot(allTreeScores)
		plt.xlabel('')
		plt.ylabel('RMSE')
		plt.title('RMSE of the trees in 100 simulation datasets')
		plt.savefig('treeScores.png')
		plt.clf()
		
		
		#4. After we have generated our malignant subclone, we can start doing real subclonal expansions.
		#For a pre-defined number of cell cycles, we ask each returned subclone to expand.
		#We keep the return value and ask this to expand again.
		#at each stage, we store the subclones that have been formed.
		#This step is repeated until we have made it through 6 cell cycles.
		
		#We can do the expansion recursively. 
	
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
			print profileEntry
			#Make a new group each time the first value is within one group
			profileValues = profileEntry.split(",")
			firstVal = profileValues[0]
			print firstVal
			if firstVal != previousVal: #append the frequency value to the current bar
				print keyVals
				print currentBar
				indices = [b + width + offset for b in np.arange(len(currentBar))]
				print indices
				plt.bar(indices, currentBar, width, alpha=opacity, label=firstVal, color=colors[i]) #this may not be the nicest color to use
				#Save the ticks where we set specific values
				[ticks.append(index) for index in indices]
				offset += len(currentBar)
				#clear the values from the previous bar group
				currentBar = []
				keyVals = []
				barWidth += 0.25
				
				print "offset: ", offset
				i += 1
			previousVal = firstVal
			currentBar.append(profile[profileEntry])
			keyVals.append(profileEntry)
		
		#Make sure that we also add the last one
		print "current bar at end: ", currentBar
		indices = [b + width + offset for b in np.arange(len(currentBar))]
		plt.bar(indices, currentBar, width, alpha=opacity, label=firstVal)
		#Save the ticks where we set specific values
		[ticks.append(index) for index in indices]
		print indices
			
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
				print eA[row][col]
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
			
			# passed = False
			# currentCycle = cycle
			# for possibleCycle in range(currentCycle, self.cellCycles):
			# 	newSubclone = subclone.expandSubclone(simulationProbabilities, self)
			# 	cycle += 1 #we do increase the cycle count so that this is not reset each time we move on
			# 	print "new cycle count: ", cycle
			# 	if newSubclone != False:
			# 		passed = True
			# 		clones.append(newSubclone)
			# 		print "viable now!"
			# 		break #we have found a good subclone!
				
			# if passed == False:
			# 	print "still unviable"
			# 	return clones #only quit if we did not find a viable subclone. 
		else:
			print "parent's name: ", subclone.name
			print "my name: ", newSubclone.name
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
		currentNoiseLevel = self.noiseLevels[0]
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
			X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
			noisedAf = X.rvs(1)[0] #take a sample
			
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
		
	def duplicateGenome(self):
		newC = []
		for c in self.C:
			oldCTumor = c.c[1]
			newC.append(C([2, oldCTumor*2]))
		newAlleles = []
		for a in range(0, len(self.A)):
			newACount = self.A[a].ACount*2
			newBCount = self.A[a].BCount*2
			#make the new names
			#the identifiers for the new alleles depends on the ACount.
			newAlleleObjects = Alleles(newACount, newBCount)
			newAlleleObjects.alleleIdentifiers += self.A[a].alleleIdentifiers
			for newAIdentifier in range(0, (self.A[a].ACount*2 - self.A[a].ACount)):
				newAlleleObjects.alleleIdentifiers.append('A' + str(uuid.uuid4()))
			for newBIdentifier in range(0, (self.A[a].BCount*2 - self.A[a].BCount)):
				newAlleleObjects.alleleIdentifiers.append('B' + str(uuid.uuid4()))	
			newAlleles.append(newAlleleObjects)
		
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
			if len(chromosomeName) < 1:
				chromosomeName = 'X'
			elif chromosomeName[0] == '12': #we skip 12
				continue 
			else:
				chromosomeName = chromosomeName[0]
			if chromosomeName not in chromosomes: #make sure that we do not add 1 for each arm but 1 for each chromosome8x
				chromosomes.append(chromosomeName)
		
		
		wholeChromosomesLost = random.sample(chromosomes, simulator.malignantWcLosses)
		#update the C and allele status of these chromosomes.
		#Each C entry corresponds to a chromosome arm. We need to link the chromosomes to the arms and lose those positions (1 copy at a time)
		#print "chromosomes lost in malignant precursor: ", wholeChromosomesLost
		previousChromosome = 0
		for chromosome in wholeChromosomesLost:
			for arm in range(0, len(possibleArms)):
				armChromosome = re.findall("(\d+)", possibleArms[arm])
				if len(armChromosome) < 1:
					if arm == (len(possibleArms)-2):
						armChromosome = 'X'
					else:
						armChromosome = 'Y'
				else:
					armChromosome = armChromosome[0]
					#print armChromosome
					#print chromosome
				if chromosome == str(armChromosome): #if this arm matches on the chromosome number
					malignantSubclone.C[arm].c[1] = malignantSubclone.C[arm].c[1] -1
					#print "identifiers: ", malignantSubclone.A[arm].alleleIdentifiers
					randomAllele = random.sample(malignantSubclone.A[arm].alleleIdentifiers, 1)[0] #select a random allele to be lost
					#select the allele from the list of identifiers. Here no somatic variants are present yet, so we can lose any copy. 
					randomAlelleIndex = malignantSubclone.A[arm].alleleIdentifiers.index(randomAllele)
					#print "random allele: ", randomAllele
					#print "match 1: ", re.match('^A.+', randomAllele)
					#print "match 2: ", re.match('^B.+', randomAllele)
					#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount. 
					if re.match('^A.+', randomAllele) is not None:
						#print "pre-del A:", malignantSubclone.A[arm].ACount
						malignantSubclone.A[arm].ACount = malignantSubclone.A[arm].ACount - 1 #we always keep allele B
						#print "pro-del A:", malignantSubclone.A[arm].ACount
						#also remove the selected allele from the list
					if re.match('^B.+', randomAllele) is not None:
						#print "pre-del B:", malignantSubclone.A[arm].ACount
						malignantSubclone.A[arm].BCount = malignantSubclone.A[arm].BCount - 1 #we always keep allele A
						#print "pro-del B:", malignantSubclone.A[arm].ACount
					del malignantSubclone.A[arm].alleleIdentifiers[randomAlelleIndex]
					
		
		#10 arm losses
		#We select which arms are lost based on their probability: there is no probability compared to gain/loss, only relative to the other positions!
		
		#Here we should not select an arm that is completely lost (only here with malignancy, otherwise the clone in the expansion is not viable).
		#Check the current C counts, if there is one with a copy number of 1 then we can no longer lose it
		filteredPossibleArms = []
		filteredNormalizedArmLossProbabilities = []
		for arm in range(0, len(possibleArms)):
			if malignantSubclone.C[arm].c[1] > (simulator.kmin): #if there is only 1 copy left, we cannot lose any further. We select for this in the malignant precursor because it always needs to be viable!
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
		
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		armsLost = np.random.choice(filteredPossibleArms, simulator.cycleArmLosses, p=filteredNormalizedArmLossProbabilities, replace=False)
		#print "arms lost in malignant subclone: ", armsLost
		#Update the lost arms in the C and A matrix
		for arm in range(0, len(armsLost)):
			#we need to find the right index, arms lost contains the values
			lostArmInd = possibleArms.index(armsLost[arm])
			malignantSubclone.C[lostArmInd].c[1] = malignantSubclone.C[lostArmInd].c[1] - 1
			
			randomAllele = random.sample(malignantSubclone.A[lostArmInd].alleleIdentifiers, 1)[0] #select a random allele to be lost
			#select the allele from the list of identifiers. No somatic variants are introduced here yet, so we select any allele to be lost. 
			randomAlelleIndex = malignantSubclone.A[lostArmInd].alleleIdentifiers.index(randomAllele)

			#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount. 
			if re.match('^A.+', randomAllele) is not None:
				malignantSubclone.A[lostArmInd].ACount = deepcopy(malignantSubclone.A[lostArmInd].ACount) - 1 #we always keep allele B
				#also remove the selected allele from the list
			if re.match('^B.+', randomAllele) is not None:
				malignantSubclone.A[lostArmInd].BCount = deepcopy(malignantSubclone.A[lostArmInd].BCount) - 1 #we always keep allele A
			del malignantSubclone.A[lostArmInd].alleleIdentifiers[randomAlelleIndex]
			
			
			
		#Fill slots for somatic variants (starting from the previous precursor, then just update the correct objects)
		#20 random SVs
		
		somaticVariantsGained = random.sample(range(0, len(malignantSubclone.somaticVariants)), simulator.malignantSVGains)
		#Choose a random allele for each somatic variant, this is after duplication of the genome so we have many copies to choose from!
		for gainedVariant in somaticVariantsGained:
			#obtain the chromosome arm that these are located on
			chromosomeArm = simulationProbabilities.somaticVariants[gainedVariant].chromosome
			#Then select a random allele from this chromosome arm
			armIndex = possibleArms.index(chromosomeArm)
			alleleACount = malignantSubclone.A[armIndex].ACount
			alleleBCount = malignantSubclone.A[armIndex].BCount
			
			randomAllele = random.sample(malignantSubclone.A[armIndex].alleleIdentifiers, 1)[0] #select a random allele to be lost
			
			#The somatic variant can be on any of these alleles.
			variantObject = malignantSubclone.somaticVariants[gainedVariant]
			variantObject.alleles.append(randomAllele)
			variantObject.value = 1
		
		#initialize an array for these variants, 1 where the somatic variant has been gained.
		#malignantSubclone.somaticVariants = np.zeros(simulator.numberOfSVs)
		#for gainedVariant in somaticVariantsGained:
		#	malignantSubclone.somaticVariants[gainedVariant] = 1
		
		#6 copies of chromosome 12p (we still need to get these values from the probability sheet)
		#Make sure that we also add allele identifiers!
		ch12pInd = possibleArms.index('12p')
		malignantSubclone.C[ch12pInd].c[1] = 6
		#We over-represent one of the alleles.
		
		randomAlleleInd = random.sample(range(0, len(malignantSubclone.A[ch12pInd].alleleIdentifiers)), 1)[0] #select a random allele to be lost, for 12p this is completely random.

		randomAllele = malignantSubclone.A[ch12pInd].alleleIdentifiers[randomAlleleInd]
		#we have 6 copies in total, starting from 2 copies of A and two copies of B. The allele that is selected gets two extra copies with this identifier.
		newAllelePrefix = 'A'
		if re.match('^A.+', randomAllele) is not None:
			malignantSubclone.A[ch12pInd].ACount = 6 - malignantSubclone.A[ch12pInd].BCount
			
			#Duplicate the random allele the same number of times
			dupAlleleCopies = (6 - malignantSubclone.A[ch12pInd].BCount)
			newAllelePrefix = 'A'
			
		if re.match('^B.+', randomAllele) is not None:
			
			malignantSubclone.A[ch12pInd].BCount = 6 - malignantSubclone.A[ch12pInd].ACount
			
			dupAlleleCopies = (6 - malignantSubclone.A[ch12pInd].ACount)
			newAllelePrefix = 'B'
		
		for dupAlleleCopy in range(0, dupAlleleCopies): #add the extra copies of the duplicated allele. Give this allele a new identifier, otherwise we get issues with somatic variants
			newAlleleName = newAllelePrefix + str(uuid.uuid4())
			malignantSubclone.A[ch12pInd].alleleIdentifiers.append(newAlleleName)

		#print "malignant precursor stats: "
		
		#for c in range(0, len(malignantSubclone.C)):
		#	print "arm: ", possibleArms[c]
		#	print "C: ", malignantSubclone.C[c].c[1]
		#	print "Acount: ", malignantSubclone.A[c].ACount
		#	print "BCount: ", malignantSubclone.A[c].BCount
		
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
		newSubclone.parent = self
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
			if newSubclone.C[arm].c[1] > (simulator.kmin-1): #we allow the loss of an arm, but if the total number of copies drops below kmin, the subclone is not viable.  
				filteredPossibleArms.append(possibleArms[arm])
				filteredNormalizedArmLossProbabilities.append(armLossProbabilities[arm])
		
		#If we remove some arms from this list, we also need to make sure that we re-normalize the probabilities. 		
		filteredNormalizedArmLossProbabilities = filteredNormalizedArmLossProbabilities / sum(filteredNormalizedArmLossProbabilities)
		
		##Determine how many somatic variants we have per allele
		allelesAndVariantCounts = dict()
		for arm in possibleArms:
			armIndex = possibleArms.index(arm)
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
					for identifier in newSubclone.A[armIndex].alleleIdentifiers:
						allelesAndVariantCounts[identifier] = []
					
		#print allelesAndVariantCounts
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
		#print "arms lost in new subclone: ", armsLost
		
		#Update the lost arms in the C and A matrix
		for arm in range(0, len(armsLost)):
			#we need to find the right index, arms lost contains the values
			lostArmInd = possibleArms.index(armsLost[arm])
			
			
			#From all alleles, check which one has the least somatic variants. If there is a tie, we choose a random one.
			leastSVAllele = ''
			currentBest = (simulator.numberOfSVs + 1) #we can never has less than the maximum number of SVs
			#print "we have identifiers: ", newSubclone.A[lostArmInd].alleleIdentifiers
			#if there is no identifier, we have no alleles left and the clone is not viable.
			if len(newSubclone.A[lostArmInd].alleleIdentifiers) < 1:
				return False
			
			newSubclone.C[lostArmInd].c[1] = newSubclone.C[lostArmInd].c[1] - 1 #only lose the copy after we are sure that we can lose the allele as well
			
			for allele in newSubclone.A[lostArmInd].alleleIdentifiers:
				#check how many somatic variants the alleles have
				if len(allelesAndVariantCounts[allele]) < currentBest:
					leastSVAllele = allele
			#randomAllele = random.sample(newSubclone.A[lostArmInd].alleleIdentifiers, 1)[0] #select a random allele to be lost
			
			#here we wish to remove the allele with the least number of somatic variants
			
			randomAlelleIndex = newSubclone.A[lostArmInd].alleleIdentifiers.index(leastSVAllele)
			
			#if the random allele matches on A, we remove the ACount. Otherwise, reduce the Bcount.
			# print "prev allele A count: ", newSubclone.A[lostArmInd].ACount
			# print "prev allele B count: ", newSubclone.A[lostArmInd].BCount
			# print "current alleles: ", newSubclone.A[lostArmInd].alleleIdentifiers
			# print "losing allele: ", leastSVAllele
			# 
			if re.match('^A.+', leastSVAllele) is not None:
				newSubclone.A[lostArmInd].ACount = newSubclone.A[lostArmInd].ACount - 1 #we always keep allele B
				#also remove the selected allele from the list
			if re.match('^B.+', leastSVAllele) is not None:
				newSubclone.A[lostArmInd].BCount = newSubclone.A[lostArmInd].BCount - 1 #we always keep allele A
			del newSubclone.A[lostArmInd].alleleIdentifiers[randomAlelleIndex]
			
			# 
			#We also remove a somatic variant if it has been lost and keep an annotation of that.
			#print "lost counts: ", allelesAndVariantCounts[leastSVAllele]
			#print "allele lost: ", leastSVAllele
			
			lostVariants = allelesAndVariantCounts[leastSVAllele]
			for lostVariant in lostVariants:
				#print lostVariant
				print "losing variant"
				lostVariant.lost = True
				lostVariant.value = 0
				
		#Gain alleles: it is always possible to gain, but only to gain alleles that have not yet been lost
		
		#If we gain the allele that contains a somatic variant, we need to duplicate these as well. 
		armsGained = np.random.choice(possibleArms, simulator.cycleArmGains, p=normalizedArmLossProbabilities, replace=False)
		#print "arms gained in new subclone: ", armsGained
		for arm in range(0, len(armsGained)):
			gainedArm = possibleArms.index(armsGained[arm])
			
			
			aAlleles = []
			bAlleles = []
			#print "allele identifiers: ", newSubclone.A[gainedArm].alleleIdentifiers
			for allele in newSubclone.A[gainedArm].alleleIdentifiers:
				#if the arm matches on A, we add it to the A alleles, otherwise the B alleles.
				
				if re.match('^A.+', allele) is not None:
					aAlleles.append(allele)
					#print "a alleles: ", aAlleles
				else:
					bAlleles.append(allele)
					#print "b alleles: ", bAlleles
			#print allele
			#print "ac: ", newSubclone.A[gainedArm].ACount
			#print "bc: ", newSubclone.A[gainedArm].BCount
			#Here we do not allow an allele to be gained when the allele has already been lost!
			newAllelePrefix = 'A'
			if newSubclone.A[gainedArm].ACount < 1 and newSubclone.A[gainedArm].BCount > 0:
				newSubclone.A[gainedArm].BCount = newSubclone.A[gainedArm].BCount + 1
				#also choose which one to gain out of the B alleles
				randomAllele = random.sample(bAlleles, 1)[0]
				#then when we have selected one, we duplicate it in the allele identifiers
				newAllelePrefix = 'B'
			elif newSubclone.A[gainedArm].BCount < 1 and newSubclone.A[gainedArm].ACount > 0:
				newSubclone.A[gainedArm].ACount = newSubclone.A[gainedArm].ACount + 1
				randomAllele = random.sample(aAlleles, 1)[0]
			elif newSubclone.A[gainedArm].BCount > 0 and newSubclone.A[gainedArm].ACount > 0: #this should generally always be the case if we make sure that the kmin is not exceeded
				#select any old allele
				#why is there sometimes no allele here? Either of the above should work. This only means that sometimes even
				#print "identifier: ", newSubclone.A[gainedArm].alleleIdentifiers
				randomAllele = random.sample(newSubclone.A[gainedArm].alleleIdentifiers, 1)[0]
				#print randomAllele
				if re.match('^A.+', randomAllele) is not None:
					newSubclone.A[gainedArm].ACount = newSubclone.A[gainedArm].ACount + 1 #we always keep allele B
				if re.match('^B.+', randomAllele) is not None:
					newSubclone.A[gainedArm].BCount = newSubclone.A[gainedArm].BCount + 1 #we always keep allele A
					newAllelePrefix = 'B'
			else: #if we have no copies left, we cannot do anything and w should skip this gain.
				continue
			newAlleleName = newAllelePrefix + str(uuid.uuid4())
			newSubclone.A[gainedArm].alleleIdentifiers.append(newAlleleName)
			#we also notify the somatic variant that it has been duplicated
			#print "random allele is: ", randomAllele
			duplicatedVariants = allelesAndVariantCounts[randomAllele]
			for duplicatedVariant in duplicatedVariants:
				duplicatedVariant.alleles.append(newAlleleName)
				
			#Only update C when we get here
			newSubclone.C[gainedArm].c[1] = newSubclone.C[gainedArm].c[1] + 1
		
		#Gain somatic variants: only gain somatic variants that have not been lost yet! We need to give these a tag to see if they have been lost or not in a previous iteration.
		#Currently losses are not implemented yet, so we can go with everything that is not already there.
		
		availableSomVarSlots = []
		for somVar in newSubclone.somaticVariants:
			#also check if the allele is still available, otherwise we cannot gain an allele there. 
			chromosome = somVar.chromosome
			chromosomeInd = possibleArms.index(somVar.chromosome)
			#print "the content: ", newSubclone.A[chromosomeInd].alleleIdentifiers
			if len(newSubclone.A[chromosomeInd].alleleIdentifiers) < 1: #if there are no alleles to add the variant to, we skip it. 
				continue
			
			if somVar.value < 1: 
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
				#it knows on which chromosome it has been gained on, but we still need to assign an allele identifier
				#if the chromosome arm already has alleles, we can use these
				#we consider the entire arm to consist of one big allele
				variantArm = somaticVariantGained.chromosome
				variantArmInd = possibleArms.index(variantArm)
				#Check the alleles that it has and select a random one
				randomAllele = random.sample(newSubclone.A[variantArmInd].alleleIdentifiers, 1)[0]
				somaticVariantGained.alleles.append(randomAllele)
			
		#Also check for viability, report a message to the method when we cannot go any further from here because the subclone is not viable (a bit of a recursive fashion).
		#Our viability checks are:
		#- We may not exceed our kmin and kmax
		#- 12p amplifications are not allowed to be lost (is this possible now?)
		i = 0
		for c in newSubclone.C:
			#print possibleArms[i]
			if simulator.chromosomeArms[i] == '12p': #exclude 12p from the viability check
				i += 1
				continue
			if c.c[1] < simulator.kmin or c.c[1] > simulator.kmax:
				return False
			i += 1
			
			
		
		#for c in range(0, len(newSubclone.C)):
		#	print "arm: ", possibleArms[c]
		#	print "C: ", newSubclone.C[c].c[1]
		#	print "Acount: ", newSubclone.A[c].ACount
		#	print "BCount: ", newSubclone.A[c].BCount
		return newSubclone
		
	
class SimulationProbabilities:
	
	somaticVariants = None #this can be a list with objects
	armProbabilities = None #make this a list of tuples to preserve sorting. 
		
	#read the file with probabilities
	
	def readFullProbabilityFile(self, file):
		self.armProbabilities = []
		self.somaticVariants = []
		#Read the file, store the somatic variant positions and arm probabilities separately
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
					self.somaticVariants.append(somVar)
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
	