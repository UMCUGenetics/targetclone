#Main script
import numpy as np
import itertools
from editDistance import EditDistance
from copy import deepcopy
import math
import random

from mu import Mu
from sample import Sample
from c import C
from fitness import Fitness
from combinations import CMuCombination
from combinations import CCombinations
from laf import LAF
from alleles import Alleles
from distance import EventDistances
from distance import PairwiseDistance
from distance import SomaticVariantDistance
from fst import FST
from sampleParser import SampleParser
from tree import Graph
from tree import SpanningArborescence
from segmentations import Segmentation

import settings #file with all tool settings
import time


#Object containing the most important starting functions of TargetClone.
class TargetClone:
	
	kmin = 1 #default values, will be overwritten by the settings
	kmax = 6
	segmentation = None
	eventDistances = None
	combinationMatrix = None
	cCombinations = None
	
	#Initialize TargetClone, this function is only called when the settings have changed and the method needs to be re-initialized. Otherwise we can load all values from history.
	#Given the kmin and kmax we can initialize the computation of the event distance between all combinations of alleles that can be made with this kmin and kmax
	#Every combination of C and mu can also be pre-initialized, including the creation of the mixture models, which is rather time consuming to do on-the-fly. 
	def __init__(self):

		self.kmin = settings.general['kmin']
		self.kmax = settings.general['kmax']
		
		#Compute event distances between all combinations of alleles
		self.eventDistances = EventDistances(self.kmin, self.kmax)
		
		#Initialize C+mu combinations and store in a matrix, objects do not need to be remade in loop
		#These are the combinations that we will explore and score for each LAF in the method
		self.combinationMatrix = np.empty([101, self.kmax], dtype='object') 
		for muIndex in range(0,101): #we allow a maximum of 0-100, so it is fine to use these settings here. 
			muClass = Mu(muIndex)
			for k in range(self.kmin, self.kmax+1):
				c = C([2, k])
				cMuCombination = CMuCombination(c, muClass, self.eventDistances)
				self.combinationMatrix[muIndex][k-1] = cMuCombination
		
		self.cCombinations = CCombinations(self.kmin, self.kmax) #make the possible combinations that can be made with C between kmin and kmax (This represents matrix Cc)
		

	#Provide the sample objects, this function will run all steps of TargetClone
	def run(self, samples):

		

		#Quick and dirty, references to self need to be updated in rest of code
		eventDistances = self.eventDistances
		cCombinations = self.cCombinations
		
		#get all somatic variants in one binary numpy array for simplicity when checking the ISA
		allSomaticVariants = self.mergeVariantsIntoMatrix(samples)
		
		#set the segmentation in each sample
		samples = self.updateSampleSegmentation(samples)
		
		measurementLength = len(samples[0].measurements.measurements)
		
		#Define the original graph, here the parent is always the healthy cell (or precursor as defined in settings) for every subclone.
		vertexNames = [] #use this to convert back from sample numbers to the actual sample names later
		for sample in samples:
			vertexNames.append(sample.name)
		vertices = range(0,len(samples))
		
		#SITUATION WHERE EACH PARENT IS 2N
		
		#Define the first graph, where node 0 (precursor) is the parent of all subclones
		edgeList = []
		for i in range(1, len(vertices)):
			edgeList.append((0, 0, i))
		currentGraph = Graph(vertices, set(edgeList), edgeList)
		
		
		
		
		#ALTERNATIVE: RANDOM PARENT FOR EACH SUBCLONE
		#Steps: - simulate random copy numbers and alleles for each parent
		# edgeList = []
		# for i in range(1, len(vertices)):
		# 	#Each parent is possible, except for itself
		# 	possibleParentsFront = range(0, i) 
		# 	possibleParentsEnd = range(i+1, len(vertices))
		# 	possibleParents = possibleParentsFront + possibleParentsEnd
		# 	
		# 	#sample a random parent
		# 	parent = random.choice(possibleParents)
		# 	edgeList.append((0, parent, i))
		# 	samples[i].parent = samples[parent]
		# currentGraph = Graph(vertices, set(edgeList), edgeList)
		# print edgeList
		
		maxIterNum = settings.general['maximumIterations'] -1
		converged = False
		iteration = 0
		#store all the graphs that we have made up until that point, this will be part of the visual output
		graphList = [] #we use a list to easily compare if a graph has been seen before in a previous iteration
		iterationGraphs = dict() #here the graphs are stored per iteration
		edges = dict()
		iterationMu = dict() #store the mu of the samples per iteration
		iterationMessages = dict()
		while converged is not True:
			
			print "Iteration ", iteration + 1
			print "Inferring the best C and mu per sample..."
			#1. Infer the best C and mu combination per sample given our current tree
			
			startTime = time.time()
			
			samples = self.inferBestCMu(samples, currentGraph)
			
			endTime = time.time()
			print "Time to infer the best C and mu for all samples: ", endTime - startTime
			
			#The best C and mu are stored for each sample, we wish to obtain this information for each sample individually and store it in a matrix that we can use to compute distances easily to infer a tree. 
			[cMatrix, aMatrix] = self.getCAndAMatrices(samples)
			
			#Compute the distance matrix (based on alleles) that we use to infer the tree
			[dists, aDists, sDists, snvAnnotations, alleleAnnotations] = self.computeDistanceMatrix(aMatrix, samples)
			
			print "Reconstructing the subclonal evolution tree"
			
			startTime = time.time()
			#Generate the initial tree
			fullGraph = self.generateInitialTree(vertexNames, dists, aDists, samples, snvAnnotations, alleleAnnotations)
			endTime = time.time()
			print "time to generate the initial tree: ", endTime - startTime
			
			
			startTime = time.time()
			#Update the tree to resolve the ISA
			[newGraph, message] = self.updateTree(fullGraph, allSomaticVariants, samples, dists)
			
			endTime = time.time()
			
			print "Time to update the tree: ", endTime - startTime
			
				
			#make a check for convergence, see if the current tree has the same edges as the newly inferred edges. Are the weights the same? Otherwise continue until the maximum allowed iterations are reached.
			#we should keep the current graphs in a list and report them sorted based on their scores.
			if iteration == maxIterNum:
				converged = True
			else:
				for previousGraph in graphList:
					if newGraph.compareIfGraphIsTheSame(previousGraph) is True:
						
						converged = True
			
			graphList.append(currentGraph)
			currentGraph = newGraph
			
			
			print "found tree: ", currentGraph.edgeList
			
			#Store the graph of each iteration
			iterationGraphs[iteration] = currentGraph
			
			#Loop through the current mu (tumor fraction) of the samples, store these in a dictionary
			#Save the mu and the messages (is ISA resolved or not) per iteration
			allMu = []
			for sample in samples:
				allMu.append(sample.bestCMu[0].mu.mu[1])
			
			iterationMu[iteration] = allMu
			iterationMessages[iteration] = message
			iteration += 1
				
			import os
			import psutil
			process = psutil.Process(os.getpid())
			print(process.memory_info().rss)
			
		#reset the node names from numbers to the actual sample names
		iterationGraphs = self.updateTreeNodeNames(iterationGraphs, vertexNames)
		
		return [cMatrix, aMatrix, samples, iterationGraphs, iterationMu, iterationMessages] #are the samples references between classes? Otherwise we may not need to return it for the mu values. 
	
	#Obtain all somatic variants from the samples and store it into one big matrix
	def mergeVariantsIntoMatrix(self, samples):
		allSomaticVariants = np.zeros((len(samples[0].somaticVariants), len(samples)))
		for sample in range(0, len(samples)):
			allSomaticVariants[:,sample] = samples[sample].somaticVariants
		return allSomaticVariants
	
	#Add the p/q segmentation to each sample
	def updateSampleSegmentation(self, samples):
		for sample in range(0, len(samples)):
			samples[sample].measurements.segmentation = self.segmentation
			
		return samples
	
	#Infer the best combination of C and mu for every measurement position in the sample with the given sample index		
	def inferBestCMuPerSample(self, samples, sampleInd, currentGraph, fitness):
		
		#we use the information in the graph to determine who is the parent of the current sample
		parentInd = currentGraph.getParentByVertexId(sampleInd)
		
		#if the parent index is larger than the number of samples, the parent is a precursor sample (if enabled in settings). THe precursor does not have allele information, so we use the parent of the precursor instead. 
		if  parentInd == len(samples): 
			parentOfPrecursor = currentGraph.getParentByVertexId(parentInd)
			samples[sampleInd].parent = samples[parentOfPrecursor]
		else:
			samples[sampleInd].parent = samples[parentInd]
		
		laf = samples[sampleInd].measurements

		bestMuScore = -float("Inf")
		samples[sampleInd].bestCMu = None #if we do not have a solution, leave it at None.
		muMin = int(settings.general['muMin'])
		muMax = int(settings.general['muMax']) + 1
		for muIndex in range(muMin,muMax):
			
			#For any other position, we need to search through 6 combinations, only for the new (i+1) position. 
			totalScore = 0
			i = 0
			bestCombinationScore = -float("Inf")
			bestCombinations = [None]*4 #placeholder for the best combinations of C that we find for this certain mu
			currentBestCMu = [None]*len(laf.measurements) #placeholder for the best C and mu combinations that we find for each measurement (same mu for each measurement)
			
			#For the first position, we need to search through 2^6 combinations (2 positions in the parent, these are known, and the same 2 in the current subclone)
			#The parent's copy numbers are fixed and known from the previous step. The child sample can have different copy numbers, but only 2^6
			#The 2^6 combinations can be made easily, combining this with the parent's information is easy to do on the fly as well.
			
			#startTime = time.time()
			
			[bestCombinations1, bestCombinations2, bestCombinationScore] = self.getBestFullCombination(samples, sampleInd, muIndex, i, fitness, bestCombinationScore)
			
			#endTime = time.time()
			#print "time to get first combination: ", endTime - startTime
			
			if bestCombinationScore is False: #there is no good combination for this mu, skip it
				continue
			#The best combinations that we find for these two measurement positions
			currentBestCMu[0] = bestCombinations1
			currentBestCMu[1] = bestCombinations2
			
			#add the score to the total score
			totalScore += bestCombinationScore
			
			#Also get the combinations for the following measurement positions, here we only need to look at 1 position at a time, 1 for the parent and 1 for the child.
			#We shift the position by one, so the previous one is always known in both the sample and the parent, and the current positions is also known in the parent,
			#thus we check only copy numbers between kmin and kmax once.
			#startTime = time.time()
			[currentBestCMu, score] = self.getBestCombinations(samples, sampleInd, muIndex, fitness, currentBestCMu)
			
			#endTime = time.time()
			#print "time to get all combinations: ", endTime - startTime
			totalScore += score
			
			if totalScore > bestMuScore: #if the score of the whole combination is better than the current best solution found, we update the C and mu combinations for all measurement positions
				samples[sampleInd].bestCMu = currentBestCMu
				bestMuScore = totalScore
	
	#Method to compute the best combination for all measurement positions after the first one. (See getBestFullCombination for the first two measurement positions).  			
	def getBestCombinations(self, samples, sampleInd, muIndex, fitness, currentBestCMu):
		#From this measurement position on we know the C and mu of the parent, and the C and mu of the previous position. So we only need to go through 6 c+mu combinations
		i = 1
		laf = samples[sampleInd].measurements
		totalScore = 0
		fst = FST()
		while i < (len(laf.measurements)-1):
			bestCombinationScore = 0
			bestCombinations = [None]*4
			
			#some laf may be NA. For i, this does not matter as it has already been inferred. i+1 can be NA. If i+1 is NA in either the parent or sample, set the C to NA.
			#Then the combinations need to be made with the previously available position. As two adjacent positions depend on each other, we keep searching for the next measurement
			#that is not NA and use this as the neighbor instead. 
			additive = 1 #the next position we will be using
			#There is an issue here if the last LAF is NA.
			breakLoop = False
			
			#we need to check the NA in the measurements and also in the previous C+mu combinations or A. We use all of these values, and if the inference of one of these failed (likely because the
			#measurement of the parent used in a previous position was NA), we cannot use this position in this iteration. 
			if np.isnan(laf.measurements[i+1]) or np.isnan(samples[sampleInd].parent.measurements.measurements[i+1]) or np.isnan(samples[sampleInd].parent.parent.measurements.measurements[i+1]) \
				or samples[sampleInd].parent.originalCMu[i] is None or samples[sampleInd].parent.originalCMu[i+1] is None \
				or samples[sampleInd].parent.parent.originalCMu[i] is None or samples[sampleInd].parent.parent.originalCMu[i+1] is None \
				or samples[sampleInd].parent.parent.A[i] is None or samples[sampleInd].parent.parent.A[i+1] is None or samples[sampleInd].parent.A[i] is None or samples[sampleInd].parent.A[i+1] is None:
				#Sometimes we cannot set a C for a sample because the laf is nan. If we then in later iterations need to re-use the C, the parent of the parent
				#will be nan too. So we need to skip this position
				#continue searching through i + x until we find a LAF which is not NA in either the sample or parent or parent of the parent
				for j in range(1, (len(laf.measurements) - i)):
					if i+j == len(laf.measurements)-1: #if we have checked all positions and we cannot find anything else, the last entries are all NA, so we do not continue.
						breakLoop = True
						break
					#The additive is the current position (i) + an additive such that we end up at a non-NA position. 
					if not (np.isnan(laf.measurements[i+j])) and not (np.isnan(samples[sampleInd].parent.measurements.measurements[i+j])) and not (np.isnan(samples[sampleInd].parent.parent.measurements.measurements[i+j])) \
					and (samples[sampleInd].parent.originalCMu[i+j] is not None) and (samples[sampleInd].parent.parent.originalCMu[i+j] is not None) \
					and samples[sampleInd].parent.parent.A[i+j] is not None and samples[sampleInd].parent.A[i+j] is not None:
						
						additive = j #so if we find a LAF that is not -1 that is at j = 5, then we add 5 positions (we start at i).
						break
					
			#If the last LAF is NA and there are no more LAF to look at, we can skip any further inference.
			if breakLoop is True:
				break
			
			#For every possible combination that we can make for 1 measurement position given the previous position and the two positions in the parent, compute a score that it will have that C between kmin and kmax.
			#The combinations of c to search through are stored in combinationMatrix. We need to look at the mixture models for the mu that we are currently investigating!
			
			previousFitnessParent = 0
			previousFitnessSample = 0
			
			#These values remain the same for each combination that we check
			
			#Sample 1
			cMuCombination1 = currentBestCMu[i]
			#Sample 2 (parent)
			cMuCombination3 = samples[sampleInd].parent.originalCMu[i]
			cMuCombination4 = samples[sampleInd].parent.originalCMu[i+additive]
			
			parentComb = str(cMuCombination3.c.c[1]) + str(cMuCombination4.c.c[1])
			
			
			#Compute the score that the LAF combination of the parent will get
			lafScore1 = fitness.computePCMuGivenLAF(cMuCombination1, laf.measurements[i], samples[sampleInd].parent.A[i])
			#For the parent we should also compute the P based on its parent
			lafScore3 = fitness.computePCMuGivenLAF(cMuCombination3, samples[sampleInd].parent.measurements.measurements[i], samples[sampleInd].parent.parent.A[i])
			lafScore4 = fitness.computePCMuGivenLAF(cMuCombination4, samples[sampleInd].parent.measurements.measurements[i+additive], samples[sampleInd].parent.parent.A[i+additive])
			
			
			pAlleles1 = samples[sampleInd].parent.A[i].getAllelesAsString()
			pAlleles2 = samples[sampleInd].parent.A[i+additive].getAllelesAsString()
			#Getting the alleles for this combination is the slowest part of the function. I don't know how to speed this up further
			alleles1 = cMuCombination1.getAllelesByLaf(samples[sampleInd].parent.A[i], laf.measurements[i]).getAllelesAsString()
			

			for combination in self.combinationMatrix[muIndex]: #there are 6 combinations to look through

				#Sample 1 combination
				
				cMuCombination2 = combination

				#Make a string version of the combinations to compute the distance using the FST. 
				if muIndex == 100: #when we have 100% normal cells, make sure that we do not compute the distance to a non-existent tumor cell, but instead to the normal cell.
					sampleComb = [2,2]
				else:
					sampleComb = str(cMuCombination1.c.c[1]) + str(cMuCombination2.c.c[1])
				
				startTime = time.time()
				
				#We wish to compute the fitness for our LAF and the next LAF given each combination
				
				lafScore2 = fitness.computePCMuGivenLAF(cMuCombination2, laf.measurements[i+additive], samples[sampleInd].parent.A[i+additive])
				
				#The total score is a combination of the individual scores
				totalLafScore = lafScore1 * lafScore2 * lafScore3 * lafScore4
				
				endTime = time.time()
			#	print "block 2 time: ", endTime - startTime
				
				startTime = time.time()
				
				alleles2 = cMuCombination2.getAllelesByLaf(samples[sampleInd].parent.A[i+additive], laf.measurements[i+additive]).getAllelesAsString()
				
				#Compute the FST score based on the alleles of the combinations
				fstScore = fst.computePcBasedOnAlleles(pAlleles1, pAlleles2, alleles1, alleles2, self.eventDistances)
				
				#Score the combinations (P(C,mu|LAF,T))
				combinationScore = totalLafScore / (fstScore + 1) #this score can not go to 0 unless the laf scores are 0
				
				if combinationScore > bestCombinationScore: #which combination of the 6 is the best for this measurement position? 
					bestCombinations[0] = cMuCombination1
					bestCombinations[1] = cMuCombination2
					bestCombinations[2] = cMuCombination3
					bestCombinations[3] = cMuCombination4
					bestCombinationScore = combinationScore
				endTime = time.time()
			#	print "block 3 time: ", endTime - startTime
			
			
			#Sanity check, if there is no solution for this mu, we should skip it. 
			if bestCombinations[1] is None:
				totalScore = -float("inf")
				badMu = True
				break
			
			if bestCombinationScore == 0: 
				totalScore = -float("inf")
				break #there is no solution for this mu
			
			#Store the best c mu for this laf measurement
			currentBestCMu[i+additive] = bestCombinations[1]
			
			#Take the log to compute the final score
			bestCombinationScore = np.log(bestCombinationScore)
			
			#Add to the total score. 
			totalScore += bestCombinationScore
		
			#Update the i to start at the position that is not NA in the next iteration.
			i = i + additive
			
		return [currentBestCMu, totalScore]
	
	#Infer the best combination of C for the first two measurement positions given this current mu. 
	def getBestFullCombination(self, samples, sampleInd, muIndex, i, fitness, bestCombinationScore):
		
		# print "CMu lengths: "
		# for sample in samples:
		# 	print sample.name
		# 	print len(sample.originalCMu)
		
		#print(process.get_memory_info()[0])
		
		laf = samples[sampleInd].measurements
		bestCombinations = [None]*4
		for combination in self.cCombinations.combinations:
		
			#get the right combinations by reference (-1 will only work with these self.kmin and self.kmax! This is a quick but dirty solution)
			#Sample 1
			cMuCombination1 = self.combinationMatrix[muIndex][(combination[0]-1)]
			cMuCombination2 = self.combinationMatrix[muIndex][(combination[1]-1)]
			
			#Sample 2 (parent)
			
		
			
			cMuCombination3 = samples[sampleInd].parent.originalCMu[i]
			cMuCombination4 = samples[sampleInd].parent.originalCMu[i+1]
			
			
			
			sampleComb = str(cMuCombination1.c.c[1]) + str(cMuCombination2.c.c[1])
			parentComb = str(cMuCombination3.c.c[1]) + str(cMuCombination4.c.c[1])

			#These are the C in each combination
			#We wish to compute the fitness for our LAF and the next LAF given each combination
			#We do not weigh the score if the measurement value is NA (we could here also search for the next position that is not NA, but since these are the first two positions, they determine what will happen
			#to the next positions due to the dependency. It appears better to ignore the NA position here). 
			if np.isnan(laf.measurements[i]) or samples[sampleInd].parent.A[i] is None:
				lafScore1 = 1
			else:
				lafScore1 = fitness.computePCMuGivenLAF(cMuCombination1, laf.measurements[i], samples[sampleInd].parent.A[i])
			if np.isnan(laf.measurements[i+1]) or samples[sampleInd].parent.A[i+1] is None:
				lafScore2 = 1
			else:
				lafScore2 = fitness.computePCMuGivenLAF(cMuCombination2, laf.measurements[i+1], samples[sampleInd].parent.A[i+1])
			#For the parent we should also compute the P based on its parent

			if np.isnan(samples[sampleInd].parent.measurements.measurements[i]) or samples[sampleInd].parent.parent.A[i] is None:
				lafScore3 = 1
			else:
				lafScore3 = fitness.computePCMuGivenLAF(cMuCombination3, samples[sampleInd].parent.measurements.measurements[i], samples[sampleInd].parent.parent.A[i])
			if np.isnan(samples[sampleInd].parent.measurements.measurements[i+1]) or samples[sampleInd].parent.parent.A[i+1] is None:
				lafScore4 = 1
			else:
				lafScore4 = fitness.computePCMuGivenLAF(cMuCombination4, samples[sampleInd].parent.measurements.measurements[i+1], samples[sampleInd].parent.parent.A[i+1])
			
			totalLafScore = lafScore1 * lafScore2 * lafScore3 * lafScore4
			
			#Also compute the FST-based score
			if muIndex == 100: #when we have 100% normal cells, make sure that we do not compute hte distance to a non-existent tumor cell.
				sampleComb = [2,2]
				
			
			#Compute the FST score based on the alleles 
			fstScore = 1
			if samples[sampleInd].parent.A[i] is None or samples[sampleInd].parent.A[i+1]is None:
				fstScore = 1
			else:
				fstScore = FST().computePcBasedOnAlleles(samples[sampleInd].parent.A[i].getAllelesAsString(), samples[sampleInd].parent.A[i+1].getAllelesAsString(),
											  cMuCombination1.getAllelesByLaf(samples[sampleInd].parent.A[i], laf.measurements[i]).getAllelesAsString(),
											  cMuCombination2.getAllelesByLaf(samples[sampleInd].parent.A[i+1], laf.measurements[i+1]).getAllelesAsString(), self.eventDistances)
				
			combinationScore = totalLafScore / (fstScore + 1)
			
			if combinationScore > bestCombinationScore: 
				bestCombinations[0] = cMuCombination1
				bestCombinations[1] = cMuCombination2
				bestCombinations[2] = cMuCombination3
				bestCombinations[3] = cMuCombination4
				bestCombinationScore = combinationScore
		#Store the best C for the first two positions
		#if the score is still -Inf, we have not found any good solution. If we cannot find a solution for the first two positions, searching for the next positions will not work as these depend on the first two.
		#So we skip this mu. This rarely happens!
		if bestCombinationScore == -float("Inf"):
			return [False,False,False]
		
		bestCombinationScore = np.log(bestCombinationScore)
		return bestCombinations[0], bestCombinations[1], bestCombinationScore
	
	#Handler to infer the best C and mu combination per sample. 	
	def inferBestCMu(self, samples, currentGraph):
		
		fitness = Fitness()
		threads = []
		#can we run the method for each sample in parallel? This would speed up the process by a lot. THe samples are not dependent on each other!
		for sampleInd in range(1, len(samples)): #Skip the first 'dummy' precursor sample
			startTime = time.time()
			self.inferBestCMuPerSample(samples,sampleInd,currentGraph,fitness)
			endTime = time.time()
			print "time per sample: ", endTime - startTime

			
		#After we have inferred the best C and mu for each sample, we can update the best C and mu in the samples. (we do not do this in the loop because if we update them immediately, the parent information of
		#the previous iteration gets lost!)
		for sample in range(0, len(samples)):
			
			samples[sample].originalCMu = samples[sample].bestCMu
			
			if samples[sample].bestCMu is None:
				measurementLength = len(samples[0].measurements.measurements)
				samples[sample].Mu = Mu(0) #assume 100% tumor
				#Set a default bestCMu in this case, we don't know the solution.
				samples[sample].bestCMu = [CMuCombination(C([2,2]), Mu(0), self.eventDistances)]*measurementLength
				samples[sample].originalCMu = samples[sample].bestCMu
			else:
				if samples[sample].bestCMu[0] is not None: #Update the C and mu combination in the actual sample
					samples[sample].Mu = samples[sample].bestCMu[0].mu
				else: #without a successfully inferred C and mu the sample is so complex it is most likely 100% tumor or contains a lot of subclones. 
					
					samples[sample].Mu = Mu(0) #assume 100% tumor, this is a guess. 
			
		return samples
	
	#When the best C and mu have been inferred, we can obtain matrices of C and A that we can use to compute distance matrices
	def getCAndAMatrices(self, samples):
		measurementLength = len(samples[0].measurements.measurements)
		cMatrix = np.empty([measurementLength, len(samples)], dtype=float)
		aMatrix = np.empty([measurementLength, len(samples)], dtype=object)
		
		#With this information we can also derived A
		for sample in range(0, len(samples)):
			
			#In the first round, our allele vectors in the samples will still be empty, so we append. In further iterations, we overwrite. 
			append = False
			if len(samples[sample].A) < 1:
				append = True
			laf = samples[sample].measurements	
			
			#For each measurement, check the inferred A and turn this into a matrix of A (as an object matrix)
			for lafInd in range(0, len(laf.measurements)): 
				#if the laf is nan, we should set the alleles to nan too. Same when the parent is nan, then we will also assign None to the sample.
				if np.isnan(samples[sample].measurements.measurements[lafInd]) or samples[sample].parent.A[lafInd] is None or samples[sample].bestCMu[lafInd] is None:
				
					if append == True: # iteration < 1
						samples[sample].A.append(None)
					else:
						samples[sample].A[lafInd] = None
				elif np.isnan(samples[sample].parent.parent.measurements.measurements[lafInd]): #for the case that the parent also has no A. 
					if append == True: # iteration < 1
						samples[sample].A.append(None)
					else:
						samples[sample].A[lafInd] = None
				else:
					parentAlleles = samples[sample].parent.A[lafInd] #get the alleles of the parent
					sampleLaf = samples[sample].measurements.measurements[lafInd] #get the alleles of the sample itself 
					
					alleles = samples[sample].bestCMu[lafInd].getAllelesByLaf(parentAlleles, sampleLaf) #get an allele object for this laf given our C+mu
					
					#if this is the first iteration, the allele matrix will be empty, so we must always append. 
					if append == True:# iteration < 1:
						samples[sample].A.append(alleles)
					else:
						samples[sample].A[lafInd] = alleles
			
			#Repeat but then for the C matrix	
			for cMu in range(0, len(samples[sample].bestCMu)):
				if samples[sample].bestCMu[cMu] is not None:
					#Check if the mu is 0% tumor fraction, then we should work with 2N as the copy number if we assume that it is normal.
					if samples[sample].bestCMu[cMu].mu.mu[1] == 0:
						cMatrix[cMu][sample] = float(samples[sample].bestCMu[cMu].c.c[0])
						aMatrix[cMu][sample] = Alleles(1,1)
					else:
						cMatrix[cMu][sample] = float(samples[sample].bestCMu[cMu].c.c[1])
						aMatrix[cMu][sample] = samples[sample].A[cMu]
				else:
					cMatrix[cMu][sample] = np.nan 
					aMatrix[cMu][sample] = np.nan
		
		return [cMatrix, aMatrix]
	
	#Compute the distance matrix based on the alleles and somatic variants. 
	def computeDistanceMatrix(self, aMatrix, samples):
		
		#Make a matrix with strings for easy event distance computation
		aCopy = deepcopy(aMatrix)
		for r in range(0, aCopy.shape[0]):
			for c in range(0, aCopy.shape[1]):
				if type(aMatrix[r,c]) is not float and aMatrix[r,c] is not None: #if this is None we cannot progress. A None means there is no allele information because the LAF is NaN. 
					aCopy[r,c] = aMatrix[r,c].alleleString
		
		#Compute the event distances using the FST
		alleleDistances = PairwiseDistance().fstAlleleDistance(aMatrix, samples)
		aDists = alleleDistances[1]
		#The FST also reports messages stating if LOH is introduced, this is annotated in the final tree
		alleleMessages = alleleDistances[0]

		#Compute the distances based on SNVs as well
		snvDistances = SomaticVariantDistance().computeSampleDistance(samples)
		sDists = snvDistances[1]
		somVarDistanceMessages = snvDistances[0] #loss/gain of SNV messages
		
		#Make a full distance matrix
		dists = np.ma.array(sDists * (aDists + 1), mask = False) #+1 to make sure that distances can never be 0, avoids later divisions by zero
		dists[np.isnan(dists)] = float("inf") #nan indicates that a relation is not possible, so we use inf. 
		
		return [dists, aDists, sDists, somVarDistanceMessages, alleleMessages]
	
	#Generate the intial tree with all edges between all subclones and their weights. We will infer a MST from this full tree. 
	def generateInitialTree(self, vertexNames, dists, aDists, samples, snvAnnotations, alleleAnnotations):
		#Generate a full tree with all the edges
		fullGraph = Graph(vertexNames, None, None)
		fullGraph.setEdgesFromDistances(dists)
		#Combine the edge annotations from the SNV and allele messages. 
		allMessages = deepcopy(snvAnnotations)
		for k in alleleAnnotations.keys():
			if k in allMessages:
				allMessages[k] = allMessages[k] + alleleAnnotations[k]
			else:
				allMessages[k] = alleleAnnotations[k]
		
		fullGraph.edgeAnnotations = allMessages

		#Remove the edges for which we are based on the allele distances confident that these are not possible (inf distance)
		for sample1Ind in range(0, len(samples)):
			for sample2Ind in range(0, len(samples)):
				if aDists[sample1Ind, sample2Ind] == float("inf"):
					fullGraph.removeEdge((0,sample1Ind, sample2Ind))
					
		return fullGraph
	
	#Update the full tree with all edges. We first infer the MST, and then start removing and replacing edges until the ISA has been resolved.
	#Every time the MST is made, we first check which edges violate the ISA. We find the violating edge with the highest distance and then remove it from the full tree.
	#From the full tree, we then again re-infer the MST and continue until the ISA has been resolved. If the ISA could not be resolved after all edges have been removed,
	#the MST with the fewest ISA violations is reported as the best solution.
	
	
	
	#What we need to do now to fix the issue of local maxima:
	#Make an additional loop, where we go through x iterations and each time remove random branches. The trees that we end up with are stored and have a score
	#associated to them. We rank the trees, and report the best one (ISA and lowest overall distance).
	#This should help ensure that we always at least find a tree for which the ISA is resolved, in cases where we now end up in a local minimum. 
	def updateTree(self, fullGraph, allSomaticVariants, samples, dists):
		

		#Copy the full graph, we iteratively update the full graph, but if there is no solution we can get back the original full graph and somatic variants (in case of the addition of precursor information)
		originalGraph = deepcopy(fullGraph) 
		originalSomaticVariants = deepcopy(allSomaticVariants)
		newSomaticVariants = deepcopy(allSomaticVariants)
		vertexNames = fullGraph.vertices
		
		if settings.trees['snvsEnabled'] == False: #If we wish to have a tree without using SNVs, do the following
			newEdges = SpanningArborescence().computeMinimumSpanningArborescence(fullGraph, newSomaticVariants) #First infer the MST
			#Create the graph object
			newVertices = deepcopy(vertexNames)
			newVertices.append(len(samples)) #add a new precursor state
			newGraph = Graph(newVertices, None, None)
			newGraph.setEdges(newEdges)
			
			newGraph.edgeAnnotations = fullGraph.edgeAnnotations
			
			#Return this tree without further updating.
			message = ''
			
			print "made graph: ", newGraph.edgeList
			
			
			return [newGraph, message]
		
		
		#Number of iterations to make in which we each time select a random branch to remove
		randIterations = settings.general['treeUpdatingIterations']

		allGraphsUnresolved = dict() #keep a list in which we store the ISA resolved and unresolved graph
		allGraphsResolved = dict()
		for iteration in range(0, randIterations):
			
			fullGraph = deepcopy(originalGraph) 
			#originalSomaticVariants = deepcopy(allSomaticVariants)
			newSomaticVariants = deepcopy(allSomaticVariants)
			#vertexNames = fullGraph.vertices
			
			#We go through the edges involved in causing ISA violation.
			#For either of these edges it may be efficient to place it elsewhere. This is the edge with the largest distance (meaning the largest difference in SNVs and edge weight in case of a tie).
			#Remove this edge from the tree
			#Re-run Edmond's algorithm but then without this edge until the ISA is resolved.
			
			unresolved = False
			resolved = False
			removedEdge = None
			seenPrecursorSamples = []
			iter = 0
			
			#1. Make a dict that will store all of the made trees when removing random branches
			#2. For each iteration, make a tree with Edmonds' and update it by removing random branches until we resolved the ISA.
			#3. After all iterations, look at the list of trees and see which one has the lowest distance and the ISA resolved. Priority: first ISA, then distance.
			#4. This is the tree that should be reported to the user with the appropiate messages attached. 
	
			#Rather than only storing the full score of the trees, we should keep a score indicating how many edges violate the ISA.
			#If we cannot solve the ISA for the tree, we report the tree which had the fewest number of violations.
			#For every list of trees associated with the # of violations, we can sort the trees by their weights. The tree with the smallest weight will be on top.
			treesAndIsaViolations = dict() #store the trees by the number of violations. We store the tree objects, these have weights associated.
			newGraph = Graph(vertexNames, set(), [])
			
			while resolved is False: #keep searching until the ISA has been resolved
				newEdges = SpanningArborescence().computeMinimumSpanningArborescence(fullGraph, newSomaticVariants) #infer the MST using edmonds' algorithm
				#maybe we should check if all edges are present. If not, then we also introduce a precursor (only if this is enabled in the settings)
				if newEdges is not False:
					precursorNeeded = False
					childrenFound = []
					for edge in newEdges:
						child = edge[2]
						for sampleInd in range(0, len(samples)):
							if child == sampleInd:
								childrenFound.append(child)
					if len(childrenFound) != (len(samples)-1): #we will always miss sample 0, this is a parent and not a child.
						#print "WARNING: the tree misses nodes, possibly an unsampled precursor is required? This can be enabled in the settings."
						precursorNeeded = True
				
				if newEdges is False and settings.trees['precursor'] is False: #if we cannot resolve a valid tree and no precursor should be added, the tree is unresolved and we report the tree with the fewest violations
					treesAndIsaViolations[float("inf")] = []
					treesAndIsaViolations[float("inf")].append(deepcopy(newGraph))
					unresolved = True
					break
				
				#We only attempt to introduce a precursor if this is specified in the settings. 
				if newEdges is False or precursorNeeded is True and settings.trees['precursor'] is True: #we end up here when there is no more possible tree. In this case, we need to reset the tree and add precursor nodes.
					[newEdges, newSomaticVariants, seenPrecursorSamples, fullGraph] = fullGraph.addPrecursorNode(originalGraph, originalSomaticVariants, samples, seenPrecursorSamples, dists)
					vertexNames.append(len(samples)) #Do this here so that we do not have to pass too many variables around
					
					if newEdges is None and newSomaticVariants is None and seenPrecursorSamples is None: #we failed to infer the best tree
						treesAndIsaViolations[float("inf")] = []
						treesAndIsaViolations[float("inf")].append(deepcopy(newGraph))
						unresolved = True
						break
				
					
				#Re-make a new tree with the newly inferred MST edges
				newVertices = deepcopy(vertexNames)
				newVertices.append(len(samples)) #add a new precursor state
				newGraph = Graph(newVertices, None, None)
				newGraph.setEdges(newEdges)
				
				newGraph.edgeAnnotations = fullGraph.edgeAnnotations
				
				#Check if the MST does not violate the ISA. If it does, bad edges are reported
				badEdgeData = newGraph.checkIfTreeIsValid(newSomaticVariants)
				badEdges = badEdgeData[0]
				
				violatingEdges = badEdgeData[1]
				
				#There are no bad edges violating the ISA left, so our tree is finished. 
				if badEdges is None:
					
					resolved = True
					#Check if the tree is not the same as in the previous iteration before we start appending
					strEdges = str(newGraph.edgeList)
					
					#Check if we have all edges. If not, then this one is actually unresolved.
					if len(childrenFound) / (len(samples)-1) < 0.8:
						allGraphsUnresolved[strEdges] = newGraph
						break
					if strEdges not in allGraphsResolved:
						#allGraphsResolved.append(newGraph)
						
						allGraphsResolved[strEdges] = newGraph
					
					break
				
				
				#The number of violating edges is a score for the trees.
				if len(violatingEdges) not in treesAndIsaViolations.keys():
					treesAndIsaViolations[len(violatingEdges)] = []
				
				#Store the trees and the number of violations. We can later report the best tree which is the one with the fewest bad edges
				treesAndIsaViolations[len(violatingEdges)].append(deepcopy(newGraph))
				#remove the problematic edge from the full tree (randomly selected)
				currentWorstEdge = random.choice(badEdges)
				
				fullGraph.removeEdge(currentWorstEdge)
			
			#Report on the best tree of this iteration.
			
			if unresolved is True or len(childrenFound) / (len(samples)-1) < 0.8: #if we did not succeed with introducing a precursor we should also report the best tree
				#print "Did not resolve the ISA, selecting the tree with the fewest violations"
				#message = "Did not resolve the ISA, reporting the tree with the fewest ISA violations"
				#in this case we select the best tree.
				bestKey = float("inf")
				
				for violationKey in treesAndIsaViolations.keys():
					if violationKey < bestKey:
						bestKey = violationKey
				
				#If the bestKey is infinite, we were unable to reconstruct a tree.
				if bestKey == float("inf"):
					bestTree = newGraph
					return bestTree
				
				#Obtain the set of trees with the least number of violations, this is the best tree for this iteration
				bestTree = None
				bestTreeWeight = float("inf")
				for tree in treesAndIsaViolations[bestKey]:
					if tree.getTotalWeight() < bestTreeWeight:
						bestTree = tree
						bestTreeWeight = tree.getTotalWeight()
				
				newGraph = deepcopy(bestTree) #use the first made graph, this is without edit operations and precursors that failed.
				strEdges = str(newGraph.edgeList)
				allGraphsUnresolved[strEdges] = newGraph

		### Checking which tree is the best, report this one
		#For each of the reconstructed trees, we first select the best (lowest distance) with the most nodes.
		#What about the ones with not all nodes reconstructed? -> leave this in unresolved
		allGraphsResolvedList = []
		allGraphsResolvedDistances = []
		reportResolved = False
		message = ""
		if len(allGraphsResolved) > 0:
			reportResolved = True
			for graph in allGraphsResolved:
				
				#calculate the distance
				treeTotalDistance = allGraphsResolved[graph].getTotalWeight()
				allGraphsResolvedDistances.append(treeTotalDistance)
				allGraphsResolvedList.append(allGraphsResolved[graph])
		
		#sort the distances, and sort the graphs accordingly
		resolvedSortedInd = np.argsort(np.array(allGraphsResolvedDistances))
		allGraphsResolvedSorted = list(np.array(allGraphsResolvedList)[resolvedSortedInd])
		
		#If we already found a good resolved one, prioritize this and report it. Otherwise check through the unresolved trees. 
		if reportResolved == True:
			print "reporting a resolved tree "
			return[allGraphsResolvedSorted[0], message]
		else: #Then repeat for the unresolved trees
			
			treesLackingNodes = []
			treesLackingNodesDistances = []
			fullGraphsUnresolvedList = []
			fullGraphsUnresolvedDistances = []
			#1. Check if there are any trees with a lower number of nodes. Push these to the back automatically.
			for graph in allGraphsUnresolved:
				
				minimumTreeContent = 0.8
				if len(childrenFound) / (len(samples)-1) < minimumTreeContent:
					#In this case, push the trees to the back.
					treesLackingNodes.append(allGraphsUnresolved[graph])
					treesLackingNodesDistances.append(allGraphsUnresolved[graph].getTotalWeight())
				else:
					fullGraphsUnresolvedList.append(allGraphsUnresolved[graph])
					fullGraphsUnresolvedDistances.append(allGraphsUnresolved[graph].getTotalWeight())
				
				
			#Sort the trees, full trees first, then the ones lacking nodes, and append.
			unresolvedFullSortedInd = np.argsort(np.array(fullGraphsUnresolvedDistances))
			fullGraphsUnresolvedSorted = np.array(fullGraphsUnresolvedList)[unresolvedFullSortedInd]
			
			#Sort the missing nodes trees
			unresolvedMissingSortedInd = np.argsort(np.array(treesLackingNodesDistances))
			lackingGraphsUnresolvedSorted = np.array(treesLackingNodes)[unresolvedMissingSortedInd]
			
			allUnresolvedGraphs = list(fullGraphsUnresolvedSorted) + list(lackingGraphsUnresolvedSorted)
			
			message = "Did not resolve the ISA, reporting the tree with the fewest ISA violations"
			if len(fullGraphsUnresolvedSorted) < 1:
				print "WARNING: tree misses nodes. Possibly an unsampled precursor is required? This can be enabled in the settings."
				print "Less than 80% of nodes are placed in the evolutionary tree. Reporting the minimum spanning tree with fewest violations instead, ISA is not resolved"
			
			return [allUnresolvedGraphs[0], message]
		# #check if the new graph contains all nodes, show the warning the the user. This will also be shown in the visualized tree output. 
		# childrenFound = []
		# for edge in newGraph.edges:
		# 	child = edge[2]
		# 	for sampleInd in range(0, len(samples)):
		# 		if child == sampleInd:
		# 			childrenFound.append(child)
		# if len(childrenFound) != (len(samples)-1): #we will always miss sample 0, this is a parent and not a child. 
		# 	print "Warning: missing too many samples to resolve the ISA, reporting truncated trees"
		
		#we also need to check if all edges are there. A tree that is truncated at only one or two positions is not bad, we can still report this to the user and not place one or two nodes.
		#if many more nodes are missing, the tree does not make much sense anymore. Here we can check if most of the tree is missing. If we miss more than 80% of nodes, report the
		#minimum spanning tree instead of an empty/half-empty tree.
		#This step does not work if the alleles are also involved in the lack of nodes! Then we need to work with a completely different set of weights. 
		
		
		return [newGraph, message]
	
	#Update the node names of the tree back to the actual sample names rather than the numbers that we use for easy edge comparison. 
	def updateTreeNodeNames(self, iterationGraphs, vertexNames):
		#We update hte node names and then save the updated tree in a dictionary with all the trees per iteration. These will later be reported to the user. 
		
		for graphInd in range(0,len(iterationGraphs)):

			graph = iterationGraphs[graphInd]
			
			#update the edge names to the sample names in the original file
			newEdgeNames = []
			for edge in iterationGraphs[graphInd].edgeList:
				
				newEdgeP = vertexNames[edge[1]]
				newEdgeC = vertexNames[edge[2]]
				newEdgeNames.append((edge[0], newEdgeP, newEdgeC))
			
			iterationGraphs[graphInd].setEdges(set(newEdgeNames))
		return iterationGraphs