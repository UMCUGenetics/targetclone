#Main script
import numpy as np
from scipy.spatial.distance import pdist, squareform
import itertools
from editDistance import EditDistance
from multiprocessing import Pool
from operator import mul
from functools import partial
import scipy.spatial.distance
#import matplotlib.pyplot as plt
from copy import deepcopy
import math

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
from tree import Kruskal
from tree import SpanningArborescence
from segmentations import Segmentation

import settings #file with all tool settings

import time
import sys

#make obj and keep things initialized
class TargetClone:
	
	kmin = 1 #default values
	kmax = 6
	segmentation = None
	eventDistances = None
	combinationMatrix = None
	cCombinations = None
	
	def __init__(self):
		self.kmin = settings.general['kmin']
		self.kmax = settings.general['kmax']
		
		#Make a dummy bestCMu for the healthy sample
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
		
		self.cCombinations = CCombinations(self.kmin, self.kmax) #re-make this, we use iterators, otherwise the combinations need to be stored in memory. 
		
	
	#Provide the sample objects, this function will run TargetClone (needs to be cleaned up and put into better objects and functions)
	def run(self, samples):

		eventDistances = self.eventDistances
		cCombinations = self.cCombinations
		
		#get all somatic variants in one binary numpy array
		allSomaticVariants = self.mergeVariantsIntoMatrix(samples)
		
		#set the segmentation
		samples = self.updateSampleSegmentation(samples)

		measurementLength = len(samples[0].measurements.measurements)
		
		#Define the original graph, here the parent is always the healthy cell for every subclone.
		vertexNames = [] #use this to convert back to the actual names later
		for sample in samples:
			vertexNames.append(sample.name)
		vertices = range(0,len(samples))
		
		#Define the first graph, where node 0 is the parent of all subclones
		edgeList = []
		for i in range(1, len(vertices)):
			edgeList.append((0, 0, i))
		currentGraph = Graph(vertices, set(edgeList), edgeList)
		
		maxIterNum = settings.general['maximumIterations'] -1
		converged = False
		iteration = 0
		graphList = [] #store all the graphs that we have made up until that point
		#Keep a table that we print to a file afterwards in which we store the edges that have been present and their occurrence at each iteration
		
		iterationGraphs = dict()
		edges = dict()
		iterationMu = dict() #store the mu of the samples per iteration
		iterationMessages = dict()
		while converged is not True:
		
			start_time = time.time()
			#This is where we get back after one iteration. The samples need to have their parents updated based on the current graph
			
			print "iteration: ", iteration
			
			#1. Infer the best C and mu combination per sample given our current tree
			samples = self.inferBestCMu(samples, currentGraph)
									
			#The best C and mu are stored for each sample, we wish to obtain this information for each sample individually and store it in a matrix that we can use to compute distances easily to infer a tree. 
			[cMatrix, aMatrix] = self.getCAndAMatrices(samples)
			
			#Compute the distance matrix that we use to infer the tree
			[dists, aDists, sDists, snvAnnotations, alleleAnnotations] = self.computeDistanceMatrix(aMatrix, samples)
			
			#Generate the initial tree
			fullGraph = self.generateInitialTree(vertexNames, dists, aDists, samples, snvAnnotations, alleleAnnotations)
			
			#Update the tree to resolve the ISA
			[newGraph, message] = self.updateTree(fullGraph, allSomaticVariants, samples)
			
				
			#make a check for convergence, see if the current tree has the same edges as the newly inferred edges. Are the weights the same? Otherwise continue until a max.
			#we should keep the current graphs in a list and report them based on their scores.
			if iteration == maxIterNum:
				converged = True
			else:
				for previousGraph in graphList:
					if newGraph.compareIfGraphIsTheSame(previousGraph) is True:
						
						converged = True
			
			graphList.append(currentGraph)
			currentGraph = newGraph
			
			# newEdgeList = []
			# for edge in newGraph.edgeList: #append the edges to the unique list, we will print this later
			# 	newEdge = (0, edge[1], edge[2])
			# 	newEdgeList.append(newEdge)
			# 	if newEdge not in uniqueEdges:
			# 		uniqueEdges.append(newEdge)
			# 	
			# edges[iteration] = newEdgeList
			
			print "best graph: "
			print currentGraph.getGraph()
			
			#print the total distance in the graph as well
			print currentGraph.getTotalWeight()
			
			#Store the graph of each iteration by its weight.
			iterationGraphs[iteration] = currentGraph
			
			#Loop through the current mu (tumor fraction) of the samples, store these in a dictionary
			allMu = []
			for sample in samples:
				allMu.append(sample.bestCMu[0].mu.mu[1])
			
			iterationMu[iteration] = allMu
			iterationMessages[iteration] = message
			iteration += 1

		iterationGraphs = self.updateTreeNodeNames(iterationGraphs, vertexNames)
		
		return [cMatrix, aMatrix, samples, iterationGraphs, iterationMu, iterationMessages] #are the samples references between classes? Otherwise we may not need to return it for the mu values. 
	
	def mergeVariantsIntoMatrix(self, samples):
		allSomaticVariants = np.zeros((len(samples[0].somaticVariants), len(samples)))
		for sample in range(0, len(samples)):
			allSomaticVariants[:,sample] = samples[sample].somaticVariants
			
			
		return allSomaticVariants
	
	def updateSampleSegmentation(self, samples):
		#also set the segmentation while we're at it
		for sample in range(0, len(samples)):
			samples[sample].measurements.segmentation = self.segmentation
			
		return samples
				
	def inferBestCMuPerSample(self, samples, sampleInd, currentGraph, fitness):
		parentInd = currentGraph.getParentByVertexId(sampleInd)
		#print "parent of ", sampleInd, " is ", parentInd
		print "sample: ", samples[sampleInd].name
		if  parentInd == len(samples): #this is the extra precursor state, currently we only support 1 extra precursor
			#get the parent of the precursor
			parentOfPrecursor = currentGraph.getParentByVertexId(parentInd)
			samples[sampleInd].parent = samples[parentOfPrecursor]
		else:
			samples[sampleInd].parent = samples[parentInd] ##is this not working correctly with the new tree?
		#print "sample ind: %d " % sampleInd
		#print "parent ind: %d" % currentGraph.getParentByVertexId(sampleInd)
		
		laf = samples[sampleInd].measurements

		bestMuScore = -float("Inf")
		samples[sampleInd].bestCMu = None #if we do not have a solution, leave it at None.
		muMin = int(settings.general['muMin'])
		muMax = int(settings.general['muMax']) + 1
		for muIndex in range(muMin,muMax):
			
			#For the first position, we need to search through 2^6 combinations
			#The parent's copy numbers are fixed and known from the previous step. The child sample can have different copy numbers, but only 2^6
			#The 2^6 combinations can be made easily, combining this with the parent's information is easy to do on the fly as well. 
			#For any other position, we need to search through 6 combinations, only for the new (i+1) position. 
			totalScore = 0
			i = 0
			bestCombinationScore = -float("Inf")
			bestCombinations = [None]*4
			currentBestCMu = [None]*len(laf.measurements)
			
			#Get the full combination for the first 4 positions
			[bestCombinations1, bestCombinations2, bestCombinationScore] = self.getBestFullCombination(samples, sampleInd, muIndex, i, fitness, bestCombinationScore)
			if bestCombinationScore is False: #there is no good combination for this mu, skip it
				continue
			currentBestCMu[0] = bestCombinations1
			currentBestCMu[1] = bestCombinations2
			
			totalScore += bestCombinationScore
			
			#Also get the combinations for the following positions, here we only need to look at 1 position at a time, 1 for the parent and 1 for the child
			[currentBestCMu, score] = self.getBestCombinations(samples, sampleInd, muIndex, fitness, currentBestCMu)
			totalScore += score
			
			if totalScore > bestMuScore: #if the score of the whole combination is better than the current best
			
				samples[sampleInd].bestCMu = currentBestCMu
				bestMuScore = totalScore
				
	def getBestCombinations(self, samples, sampleInd, muIndex, fitness, currentBestCMu):
		#From here we know the parent, and the previous position. So we only need to go through 6 c mu combinations
		i = 1
		laf = samples[sampleInd].measurements
		totalScore = 0
		while i < (len(laf.measurements)-1):
		
			bestCombinationScore = 0
			bestCombinations = [None]*4
			
			#some laf may be NA. For i, this does not matter as it has already been inferred. i+1 can be NA. If i+1 is NA in either the parent or sample, set the C to NA.
			#Then the combinations need to be made with the previously available position. If the previous position is on a different chromosome, we actually need to go through
			#the 2^6 combinations step again! Then continue. For now to get things started, we can consider the genome as a long string. 
			additive = 1 #the next position we will be using
			#There is an issue here if the last LAF is NA.
			breakLoop = False
		
			if np.isnan(laf.measurements[i+1]) or np.isnan(samples[sampleInd].parent.measurements.measurements[i+1]) or np.isnan(samples[sampleInd].parent.parent.measurements.measurements[i+1]) \
				or samples[sampleInd].parent.originalCMu[i] is None or samples[sampleInd].parent.originalCMu[i+1] is None \
				or samples[sampleInd].parent.parent.originalCMu[i] is None or samples[sampleInd].parent.parent.originalCMu[i+1] is None \
				or samples[sampleInd].parent.parent.A[i] is None or samples[sampleInd].parent.parent.A[i+1] is None or samples[sampleInd].parent.A[i] is None or samples[sampleInd].parent.A[i+1] is None:
				#Sometimes we cannot set a C for a sample because the laf is nan. If we then in later iterations need to re-use the C, the parent of the parent
				#will be nan too. So we need to skip this position (for now, maybe later we can try to find a full new combination based on the laf)

				#continue searching through i + x until we find a LAF which is not -1 in either the sample or parent or parent of the parent
				for j in range(1, (len(laf.measurements) - i)):
					if i+j == len(laf.measurements)-1: #if we have checked all positions and we cannot find anything else, the last entries are all NA.
						breakLoop = True
						break
					
					if not (np.isnan(laf.measurements[i+j])) and not (np.isnan(samples[sampleInd].parent.measurements.measurements[i+j])) and not (np.isnan(samples[sampleInd].parent.parent.measurements.measurements[i+j])) \
					and (samples[sampleInd].parent.originalCMu[i+j] is not None) and (samples[sampleInd].parent.parent.originalCMu[i+j] is not None) \
					and samples[sampleInd].parent.parent.A[i+j] is not None and samples[sampleInd].parent.A[i+j] is not None:
						
						additive = j #so if we find a LAF that is not -1 that is at j = 5, then we add 5 positions (we start at i).
						break
					#We end here if the parent or sample has a NA as LAF. We cannot infer a C here, so we should set this C to NA as well.
					
				#If the last LAF is NA and there are no more LAF to look at, we can skip this step.
			if breakLoop is True:
				break

			for combination in self.combinationMatrix[muIndex]: #there are 6 combinations to look through
				#Sample 1
				cMuCombination1 = currentBestCMu[i] #Sometimes ther is no i in here
				cMuCombination2 = combination

				#Sample 2 (parent)
				cMuCombination3 = samples[sampleInd].parent.originalCMu[i]
				cMuCombination4 = samples[sampleInd].parent.originalCMu[i+additive]
				
				#print samples[sampleInd].parent.A.alleleString
				#print samples[sampleInd].parent.parent.A.alleleString

				sampleComb = str(cMuCombination1.c.c[1]) + str(cMuCombination2.c.c[1])
				parentComb = str(cMuCombination3.c.c[1]) + str(cMuCombination4.c.c[1])
				if muIndex == 100: #when we have 100% normal cells, make sure that we do not compute hte distance to a non-existent tumor cell.
					sampleComb = [2,2]
					
				# if sampleInd == 2: #we are intersted in GCNIS
				# 	print "parent combination: ", parentComb
				# 	print "sample combination: ", sampleComb
				# 	print "laf1: ", laf.measurements[i]
				# 	print "laf2: ", laf.measurements[i+additive]
				# 	print "laf3: ", samples[sampleInd].parent.measurements.measurements[i]
				# 	print "laf4: ", samples[sampleInd].parent.measurements.measurements[i+additive]
				# 	
				#These are the C in each combination
				#We wish to compute the fitness for our LAF and the next LAF given each combination
				# print "laf: ", laf.measurements[i]
				# print "combination: ", cMuCombination1.c.c[1]
				lafScore1 = fitness.computePCMuGivenLAF(cMuCombination1, laf.measurements[i], samples[sampleInd].parent.A[i])
				# print "laf 2: ", laf.measurements[i+additive]
				# print "combination 2: ", cMuCombination2.c.c[1]
				lafScore2 = fitness.computePCMuGivenLAF(cMuCombination2, laf.measurements[i+additive], samples[sampleInd].parent.A[i+additive])
				#For the parent we should also compute the P based on its parent
				# print "laf 3: ", samples[sampleInd].parent.measurements.measurements[i]
				# print "combination 3: ", cMuCombination3.c.c[1]
				lafScore3 = fitness.computePCMuGivenLAF(cMuCombination3, samples[sampleInd].parent.measurements.measurements[i], samples[sampleInd].parent.parent.A[i])
				# print "laf 4: ", samples[sampleInd].parent.measurements.measurements[i+additive]
				# print "combination 4: ", cMuCombination4.c.c[1]
				lafScore4 = fitness.computePCMuGivenLAF(cMuCombination4, samples[sampleInd].parent.measurements.measurements[i+additive], samples[sampleInd].parent.parent.A[i+additive])
				#if samples[sampleInd].parent.parent.A[i].alleleString is not 'AB':
					#print "current parental allele: %s" % samples[sampleInd].parent.parent.A[i].alleleString
				#print "next parental allele: %s" % samples[sampleInd].parent.parent.A[i+1].alleleString
				# 
				# print "combination 1: ", cMuCombination1.c.c[1]
				# print "combination 2: ", cMuCombination2.c.c[1]
				# print "combination 3: ", cMuCombination3.c.c[1]
				# print "combination 4: ", cMuCombination4.c.c[1]
				# 
				# print "Laf 1: ", laf.measurements[i]
				# print "laf 2: ", laf.measurements[i+additive]
				# print "laf 3: ", samples[sampleInd].parent.measurements.measurements[i]
				# print "laf 4: ", samples[sampleInd].parent.measurements.measurements[i+additive]
				# 
				# 
				# 

				totalLafScore = lafScore1 * lafScore2 * lafScore3 * lafScore4
				# print "total laf score: ", totalLafScore

				#Also compute the FST-based score
				fstScore = FST().computeDistance(sampleComb, parentComb)
				#what are the alleles of the combinations?
				# print "alleles parent, 1: ", samples[sampleInd].parent.A[i].getAllelesAsString()
				# print "alleles parent, 2: ", samples[sampleInd].parent.A[i+additive].getAllelesAsString()
				# print "alleles 1: ", cMuCombination1.getAllelesByLaf(samples[sampleInd].parent.A[i], laf.measurements[i]).getAllelesAsString()
				# print "alleles 2: ", cMuCombination2.getAllelesByLaf(samples[sampleInd].parent.A[i+additive], laf.measurements[i+additive]).getAllelesAsString()
				# #These are the alleles that we would choose based on this combination
				fstScore = FST().computePcBasedOnAlleles(samples[sampleInd].parent.A[i].getAllelesAsString(), samples[sampleInd].parent.A[i+additive].getAllelesAsString(),
											  cMuCombination1.getAllelesByLaf(samples[sampleInd].parent.A[i], laf.measurements[i]).getAllelesAsString(),
											  cMuCombination2.getAllelesByLaf(samples[sampleInd].parent.A[i+additive], laf.measurements[i+additive]).getAllelesAsString())
				
				
				#print "fst score: ", fstScore
				combinationScore = totalLafScore / (fstScore + 1) #this score can not go to 0 unless the laf scores are 0
				# print "combination score: ", combinationScore
				#combinationScore = totalLafScore
				# 
				
				# print "fst score: ", fstScore
				# print "total score: ", combinationScore
			
					
				if combinationScore > bestCombinationScore: #which combination of the 6 is the best? 
					bestCombinations[0] = cMuCombination1
					bestCombinations[1] = cMuCombination2
					bestCombinations[2] = cMuCombination3
					bestCombinations[3] = cMuCombination4
					bestCombinationScore = combinationScore
			
			#Sanity check, if there is no solution for this mu, we should skip it. 
			if bestCombinations[1] is None:
				totalScore = -float("inf")
				badMu = True
				break
				#continue
			
			if bestCombinationScore == 0: #we now end up here quite a lot
				print "bad combination: "
				# print bestCombinations[0].c.c[1]
				# print bestCombinations[1].c.c[1]
				# print bestCombinations[2].c.c[1]
				# print bestCombinations[3].c.c[1]
				# #Also print the LAF, what is happening?
				# print "measurements: "
				# print laf.measurements[i]
				# print laf.measurements[i+additive]
				# 
				# #print the scores, see which one is causing the score to go to infinite
				# print "score 1: ", lafScore1
				# print "score 2: ", lafScore2
				# print "score 3: ", lafScore3
				# print "score 4: ", lafScore4
				totalScore = -float("inf")
				break #there is no solution for this mu
				#scontinue #we cannot do anything with this position, this is most likely because the noise is very big. We should not quit the method here based
				#on this position alone, but let it continue searching through the rest of the positions. This allows errors in the data which are likely corrected for
				#in other positions. 
			
			#Store the best c mu for this laf 
			#currentBestCMu.append(bestCombinations[1])
			currentBestCMu[i+additive] = bestCombinations[1]
			#if bestCombinationScore == 0:
			#	print sampleInd
			#	print bestCombinations[1]
			bestCombinationScore = np.log(bestCombinationScore)
			
			totalScore += bestCombinationScore
		
			#totalScore *= bestCombinationScore
			#Update the i to start at the position that is not NA in the next iteration.
			i = i + additive
		return [currentBestCMu, totalScore]
				
	def getBestFullCombination(self, samples, sampleInd, muIndex, i, fitness, bestCombinationScore):
		laf = samples[sampleInd].measurements
		bestCombinations = [None]*4
		for combination in self.cCombinations.combinations:
			#print combination
			#get the right combinations by reference (-1 will only work with these self.kmin and self.kmax! This is a quick but dirty solution)
			#Sample 1
			cMuCombination1 = self.combinationMatrix[muIndex][(combination[0]-1)]
			cMuCombination2 = self.combinationMatrix[muIndex][(combination[1]-1)]
			#Sample 2 (parent): this is known and should be stored in the parent sample information. 
			#cMuCombination3 = combinationMatrix[muIndex][(combination[2]-1)]
			#cMuCombination4 = combinationMatrix[muIndex][(combination[3]-1)]
			
			cMuCombination3 = samples[sampleInd].parent.originalCMu[i]
			cMuCombination4 = samples[sampleInd].parent.originalCMu[i+1]
			
			sampleComb = str(cMuCombination1.c.c[1]) + str(cMuCombination2.c.c[1])
			parentComb = str(cMuCombination3.c.c[1]) + str(cMuCombination4.c.c[1])

			#These are the C in each combination
			#We wish to compute the fitness for our LAF and the next LAF given each combination
			#We do not weigh the score if the measurement value is NA.
	
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
			# 
			# print "laf 1: ", laf.measurements[i]
			# print "laf 2: ", laf.measurements[i+1]
			# print "laf 3: ", samples[sampleInd].parent.measurements.measurements[i]
			# print "laf 4: ", samples[sampleInd].parent.measurements.measurements[i+1]
			# print "combination 1: ", cMuCombination1.c.c[1]
			# print "combination 2: ", cMuCombination2.c.c[1]
			# print "combination 3: ", cMuCombination3.c.c[1]
			# print "combination 4: ", cMuCombination4.c.c[1]
			# 
			totalLafScore = lafScore1 * lafScore2 * lafScore3 * lafScore4
			# print "laf score: ", totalLafScore
			
			#Also compute the FST-based score
			if muIndex == 100: #when we have 100% normal cells, make sure that we do not compute hte distance to a non-existent tumor cell.
				sampleComb = [2,2]
				
			#This does not work if either i or i+1 is NA!
			fstScore = FST().computeDistance(sampleComb, parentComb)
			if samples[sampleInd].parent.A[i] is None or samples[sampleInd].parent.A[i+1]is None:
				fstScore = 1
			else:
				fstScore = FST().computePcBasedOnAlleles(samples[sampleInd].parent.A[i].getAllelesAsString(), samples[sampleInd].parent.A[i+1].getAllelesAsString(),
											  cMuCombination1.getAllelesByLaf(samples[sampleInd].parent.A[i], laf.measurements[i]).getAllelesAsString(),
											  cMuCombination2.getAllelesByLaf(samples[sampleInd].parent.A[i+1], laf.measurements[i+1]).getAllelesAsString())
				
			# print "fst score: ", fstScore
			
			combinationScore = totalLafScore / (fstScore + 1)
			
			if combinationScore > bestCombinationScore: 
				bestCombinations[0] = cMuCombination1
				bestCombinations[1] = cMuCombination2
				bestCombinations[2] = cMuCombination3
				bestCombinations[3] = cMuCombination4
				bestCombinationScore = combinationScore
		#Store the best C for the first two positions
		#Unless the score is 0, then we do not have a working solution
		#if bestCombinationScore == 0:
		#	continue
		if bestCombinationScore == -float("Inf"):
			return [False,False,False]
		#currentBestCMu.append(bestCombinations[0])
		#currentBestCMu.append(bestCombinations[1])
		
		bestCombinationScore = np.log(bestCombinationScore)
		return bestCombinations[0], bestCombinations[1], bestCombinationScore
		
	def inferBestCMu(self, samples, currentGraph):
		#import threading
		import time
		start_time = time.time()
		fitness = Fitness()
		threads = []
		#can we run the method for each sample in parallel? This would speed up the process by a lot. THe samples are not dependent on each other!
		for sampleInd in range(1, len(samples)): #Skip the first 'dummy' precursor sample
			self.inferBestCMuPerSample(samples,sampleInd,currentGraph,fitness)
			#t = threading.Thread(target=self.inferBestCMuPerSample, args=[samples,sampleInd,currentGraph, fitness])
			#threads.append(t)
			#t.start()
		#for thread in threads:
		#	thread.join()
		
		#After we have inferred the best C and mu for each sample, we can update the best C and mu in the samples. 
		for sample in range(0, len(samples)):
			
			samples[sample].originalCMu = samples[sample].bestCMu
			if samples[sample].bestCMu is None:
				measurementLength = len(samples[0].measurements.measurements)
				samples[sample].Mu = Mu(0) #assume 100% tumor
				#Set a default bestCMu in this case, we don't know the solution.
				print "sample ", sample, " setting CMu to 2"
				samples[sample].bestCMu = [CMuCombination(C([2,2]), Mu(0), self.eventDistances)]*measurementLength
			else:
				if samples[sample].bestCMu[0] is not None:
					print "setting mu to: ", samples[sample].bestCMu[0].mu.mu, " in sample: ", samples[sample].name
					samples[sample].Mu = samples[sample].bestCMu[0].mu
				else: #without a successfully inferred C and mu the sample is so complex it is most likely 100% tumor or contains a lot of subclones. 
					print sample, " Assuming 100% tumor"
					
					samples[sample].Mu = Mu(0) #assume 100% tumor

		print("--- %s seconds for all samples of 1 patient ---" % (time.time() - start_time))	
		return samples
	
	def getCAndAMatrices(self, samples):
		measurementLength = len(samples[0].measurements.measurements)
		cMatrix = np.empty([measurementLength, len(samples)], dtype=float)
		aMatrix = np.empty([measurementLength, len(samples)], dtype=object)
		
		#With the best C and mu inferred, we can turn it into a matrix
		#With this information we can also derived A
		for sample in range(0, len(samples)):
			
			append = False
			if len(samples[sample].A) < 1:
				append = True
			laf = samples[sample].measurements	
			#set alleles
			for lafInd in range(0, len(laf.measurements)): #what is the best place to do this?
				#if the laf is nan, we should set the alleles to nan too. Same when the parent is nan, then we will also assign None to the sample.
				
				if np.isnan(samples[sample].measurements.measurements[lafInd]) or samples[sample].parent.A[lafInd] is None or samples[sample].bestCMu[lafInd] is None:
				
					if append == True: # iteration < 1:
						samples[sample].A.append(None)
					else:
						samples[sample].A[lafInd] = None
				elif np.isnan(samples[sample].parent.parent.measurements.measurements[lafInd]): #for the case that the parent also has no A. 
					if append == True: # iteration < 1:
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
				
			for cMu in range(0, len(samples[sample].bestCMu)):
				if samples[sample].bestCMu[cMu] is not None:
					cMatrix[cMu][sample] = float(samples[sample].bestCMu[cMu].c.c[1])
					aMatrix[cMu][sample] = samples[sample].A[cMu]
				else:
					cMatrix[cMu][sample] = np.nan 
					aMatrix[cMu][sample] = np.nan
		
		return [cMatrix, aMatrix]
	
	def computeDistanceMatrix(self, aMatrix, samples):
		
		aCopy = deepcopy(aMatrix)
		for r in range(0, aCopy.shape[0]):
			for c in range(0, aCopy.shape[1]):
				if type(aMatrix[r,c]) is not float and aMatrix[r,c] is not None: #if this is None we cannot progress. A None means there is no allele information because the LAF is NaN. 
					aCopy[r,c] = aMatrix[r,c].alleleString
		
		
		alleleDistances = PairwiseDistance().fstAlleleDistance(aMatrix, samples)
		aDists = alleleDistances[1]
		alleleMessages = alleleDistances[0]
		
		aStringMatrix = deepcopy(aMatrix)
		aStringMatrix = aStringMatrix.astype(str)
		for row in range(0,aStringMatrix.shape[0]):
			for col in range(0, aStringMatrix.shape[1]):
				if aMatrix[row][col] is np.nan or aMatrix[row][col] is None:
					aStringMatrix[row][col] = 'NA'
				else:
					aStringMatrix[row][col] = aMatrix[row][col].alleleString
		
		print aDists
		
		snvDistances = SomaticVariantDistance().computeSampleDistance(samples)
		sDists = snvDistances[1]
		somVarDistanceMessages = snvDistances[0]
		print sDists
		dists = np.ma.array(sDists * (aDists + 1), mask = False) #+1 to make sure that distances can never be 0
		dists[np.isnan(dists)] = float("inf")
		print dists
		#dists.mask[np.diag_indices(dists.shape[0])] = True
		return [dists, aDists, sDists, somVarDistanceMessages, alleleMessages]
	
	def generateInitialTree(self, vertexNames, dists, aDists, samples, snvAnnotations, alleleAnnotations):
		#Generate a full tree with all the edges
		fullGraph = Graph(vertexNames, None, None)
		fullGraph.setEdgesFromDistances(dists)
		#Combine the edge annotations:
		allMessages = deepcopy(snvAnnotations)
		for k in alleleAnnotations.keys():
			if k in allMessages:
				allMessages[k] = allMessages[k] + alleleAnnotations[k]
			else:
				allMessages[k] = alleleAnnotations[k]
		
		fullGraph.edgeAnnotations = allMessages
		print fullGraph
		print fullGraph.edgeList
		#Remove the edges for which we are based on the allele distances confident that these are not possible
		for sample1Ind in range(0, len(samples)):
			for sample2Ind in range(0, len(samples)):
				if aDists[sample1Ind, sample2Ind] == float("inf"):
					fullGraph.removeEdge((0,sample1Ind, sample2Ind))
					
		return fullGraph
	
	def updateTree(self, fullGraph, allSomaticVariants, samples):
		#Copy the full graph, we iteratively update the full graph, but if there is no solution we can get back the original full graph and somatic variants (in case of precursors)
		originalGraph = deepcopy(fullGraph) 
		originalSomaticVariants = deepcopy(allSomaticVariants)
		newSomaticVariants = deepcopy(allSomaticVariants)
		vertexNames = fullGraph.vertices
		
		#We go through the edges involved in causing ISA violation.
		#For either of these edges it may be efficient to place it elsewhere. This is the edge with the largest distance.
		#Remove this edge from the tree
		#Re-run Edmond's algorithm but then without this edge until the ISA is resolved.
		
		unresolved = False
		resolved = False
		removedEdge = None
		seenPrecursorSamples = []
		edgeGone = False
		iter = 0
		savedGraph = None
		#print "resolving the ISA"
		print "reconstructing tree: "
		#Rather than only storing the full score of the trees, we should keep a score indicating how many edges violate the ISA.
		#If we cannot solve the ISA for the tree, we report the tree which had the fewest number of violations.
		#For every list of trees associated with the # of violations, we can sort the trees by their weights. The tree with the smallest weight will be on top.
		treesAndIsaViolations = dict() #store the trees by the number of violations. We store the tree objects, these have weights associated.
		newGraph = Graph(vertexNames, set(), [])
		
		#In the case that we ignore SNVs, we do not need to resolve the ISA by checking bad edges. The first found tree is the solution. 
		
		while resolved is False:
			#print fullGraph.getGraph()
			newEdges = SpanningArborescence().computeMinimumSpanningArborescence(fullGraph, newSomaticVariants)

			#maybe we should check if all edges are present. If not, then we also introduce a precursor
			if newEdges is not False:
				precursorNeeded = False
				childrenFound = []
				for edge in newEdges:
					child = edge[2]
					for sampleInd in range(0, len(samples)):
						if child == sampleInd:
							childrenFound.append(child)
				if len(childrenFound) != (len(samples)-1): #we will always miss sample 0, this is a parent and not a child.
					print newEdges
					print "Warning: the tree misses nodes"
					precursorNeeded = True
			
			if newEdges is False and settings.trees['precursor'] is False:
				treesAndIsaViolations[float("inf")] = []
				treesAndIsaViolations[float("inf")].append(deepcopy(newGraph))
				unresolved = True
				break
			
			#We only attempt to introduce a precursor if this is specified in the settings. 
			if newEdges is False or precursorNeeded is True and settings.trees['precursor'] is True: #we end up here when there is no more possible tree. In this case, we need to reset the tree and add precursor nodes.
				[newEdges, newSomaticVariants] = fullGraph.addPrecursorNode(originalGraph, originalSomaticVariants, samples)
				
			
			newVertices = deepcopy(vertexNames)
			newVertices.append(len(samples)) #add a new precursor state
			newGraph = Graph(newVertices, None, None)
			newGraph.setEdges(newEdges)
			newGraph.edgeAnnotations = fullGraph.edgeAnnotations
			
			if iter == 0:
				savedGraph = deepcopy(newGraph)
			iter += 1
			
			
			badEdgeData = newGraph.checkIfTreeIsValid(newSomaticVariants)
			badEdges = badEdgeData[0]
			violatingEdges = badEdgeData[1]
			
			#In this case, we do not want to resolve the ISA and the first reported tree is the minimum distance tree.
			#Thus, the bad edges are None by default.
			badEdges = None
			if badEdges is None:
				
				resolved = True
				break
			
			
			#THe number of violating edges is a score for the trees.
			
			if len(violatingEdges) not in treesAndIsaViolations.keys():
				treesAndIsaViolations[len(violatingEdges)] = []
			
			treesAndIsaViolations[len(violatingEdges)].append(deepcopy(newGraph))
			print "bad edges: ", badEdges
			#Remove the edge with the largest distance
			#we choose the edge that together (somvar * distance) has the worst score. 
			currentLargestDist = -float("inf")
			currentWorstEdge = 0
			if len(badEdges) > 0:
				edgeCounter = 0
				for edge in badEdges:
					child = edge[2]
					parent = edge[1]
					print "current edge: ", edge
					#totalDistance = (math.exp(violatingWeights[edgeCounter])) * (edge[0])
					totalDistance = edge[0]
					print "total distance: ", totalDistance
					if totalDistance > currentLargestDist:
						currentWorstEdge = edge
						currentLargestDist = totalDistance
						
					edgeCounter += 1
			#remove the problematic edge
			
			print "removing edge: ", currentWorstEdge
			fullGraph.removeEdge(currentWorstEdge)
			
			
			#print "removing edge: ", currentWorstEdge
		
		#if newGraph is None: #sometimes we cannot resolve the ISA
	#	newGraph = deepcopy(fullGraph)
		
		#check if the new graph contains all nodes, throw this warning only at the end 
		childrenFound = []
		for edge in newGraph.edges:
			child = edge[2]
			for sampleInd in range(0, len(samples)):
				if child == sampleInd:
					childrenFound.append(child)
		if len(childrenFound) != (len(samples)-1): #we will always miss sample 0, this is a parent and not a child. 
			print "Warning: missing too many samples to resolve the ISA, reporting truncated trees"
		
		#we also need to check if all edges are there. A tree that is truncated at only one or two positions is not bad, we can still report this to the user and not place one or two nodes.
		#if many more nodes are missing, the tree does not make much sense anymore. Here we can check if most of the tree is missing. If we miss more than 80% of nodes, report the
		#minimum spanning tree instead of an empty/half-empty tree.
		#This step does not work if the alleles are also involved in the lack of nodes! Then we need to work with a completely different set of weights. 
		minimumTreeContent = 0.8
		if len(childrenFound) / (len(samples)-1) < 0.8:
			print "Less than 80% of nodes are placed in the evolutionary tree. Reporting the minimum spanning tree instead, ISA is not resolved"
		message = ""
		if unresolved is True or len(childrenFound) / (len(samples)-1) < 0.8: #if we did not succeed with introducing a precursor we should also report the best tree
			print "Did not resolve the ISA, selecting the tree with the fewest violations"
			message = "Did not resolve the ISA, reporting the tree with the fewest ISA violations"
			#in this case we select the best tree.
			bestKey = float("inf")
			
			for violationKey in treesAndIsaViolations.keys():
				if violationKey < bestKey:
					bestKey = violationKey
			
			#If the bestKey is infinite, we were unable to reconstruct a tree.
			if bestKey == float("inf"):
				bestTree = newGraph
				print "the best tree: ", bestTree.edgeList
				print "Did not find a correct tree"
				return bestTree
			
			
			#Obtain the set of trees with this number of violations
			bestTree = None
			bestTreeWeight = float("inf")
			for tree in treesAndIsaViolations[bestKey]:
				if tree.getTotalWeight() < bestTreeWeight:
					bestTree = tree
					bestTreeWeight = tree.getTotalWeight()
			
			#print "number of violations: ", bestKey	
			
			newGraph = deepcopy(bestTree) #use the first made graph, this is without edit operations and precursors that failed.
		return [newGraph, message]
	
	def updateTreeNodeNames(self, iterationGraphs, vertexNames):
		#For every iteration we made a tree, report all of these to the user.
		#We need to do something more in the case that we cannot reconstruct the tree. This means that the precursor addition & edit operations are not sufficient, so
		#perhaps we need to start storing the minimum distance trees instead of the truncated trees.
		#Instead of the node indices we should also print the sample names. 
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