"""
	The goal of this script is to compute how often the distance from a subclone to its parent is smaller than the distance to the other subclones. 

"""
import sys
sys.path.insert(0, '../../TargetClone/')

import os
import re
import computeTreeErrorOtherMetrics
from tree import Graph
import numpy as np
from distance import EventDistances
from combinations import AlleleCombination
from alleles import Alleles
from mu import Mu

simulationFolder = "Results/generic_noise0.02/"
eventDistanceObj = EventDistances(0, 6)

#From the simulation data, get the alleles and the tree from a file.
def readSimulationData(dataFolder):
	
	trees = []
	alleles = []
	simulationDataIds = []
	
	for subdir, dirs, files in os.walk(simulationFolder):
		
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		simulationDataIds.append(subdir)
		for file in files:
		
			#Also collect the real tree and inferred tree to compute the anscestry swap errors
			if re.match('RealTrees', file): #read the file and obtain the error
				stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
				tree = eval(stringDict)
				realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
				trees.append(realTree)
			if re.match('RealA', file):
				currentAlleles = np.loadtxt(subdir + '/' + file, dtype='str')
				alleles.append(currentAlleles)
			
		
	return trees, alleles, simulationDataIds

#For every subclone, compute the distance to its parent and to the other subclones. 
def computeDistances(trees, alleles, simulationDataIds):
	
	smallerThanParentCounts = []
	
	perSubcloneCounts = dict() #keep per set per subclone how many times other subclone distances are larger than to the parent. 
	
	for simulationDataSet in range(0, len(trees)):
		
		perSubcloneCounts[simulationDataIds[simulationDataSet]] = dict()
		
		print simulationDataSet
		
		currentTree = trees[simulationDataSet]
		currentAlleles = alleles[simulationDataSet]
		currentSmallerThanParentCounts = []
		
		#For every subclone in the vertices list, determine who the parent is.
		smallerThanParentCount = 0
		vertexCombinations = 0
		print currentTree.vertices
		
		for vertexInd in range(0, len(currentTree.vertices)):
			smallerThanParentCountPerSubclone = 0
			
			vertex = currentTree.vertices[vertexInd]
			parent = ""
			for edge in currentTree.edgeList:
				if edge[2] == vertex:
					parent = edge[1]
			
			if parent == "": #skip the subclones that have no parent
				continue
				
			
					
			parentInd = currentTree.vertices.index(parent)
			
			#Compute the allelic event distance from the subclone to the parent
			subcloneAlleles = currentAlleles[:,vertexInd]
			parentAlleles = currentAlleles[:,parentInd]
			
			#Compute the event distance
			distanceToParent = 0
			for position in range(0, len(subcloneAlleles)):
				eventDistance = eventDistanceObj.getRawEventDistance(parentAlleles[position], subcloneAlleles[position]) #should use the settings
				
				distanceToParent += eventDistance
			print "vertex ", vertex, " to ", parent, " = ", distanceToParent
			#Now do the same for the non-parent subclones, and compute how often the distance is smaller than the distance to the parent.
			
			
			otherSubcloneCount = 0
			for vertexInd2 in range(0, len(currentTree.vertices)):
				vertex2 = currentTree.vertices[vertexInd2]
				
				if vertex2 == vertex:
					continue
			
				if vertex2 == parent:
					continue
				
				vertexCombinations += 1
				otherSubcloneCount += 1
				
				totalDistance = 0
				otherSubcloneAlleles = currentAlleles[:,vertexInd2]
				for position in range(0, len(otherSubcloneAlleles)):
					eventDistance = eventDistanceObj.getRawEventDistance(otherSubcloneAlleles[position], subcloneAlleles[position]) #should use the settings
					
					totalDistance += eventDistance
					
				print "vertex ", vertex, " to ", vertex2, " = ", totalDistance	
				if totalDistance < distanceToParent:
					print "vertex ", vertex, " is more similar to ", vertex2, " than to ", currentTree.vertices[parentInd]
					smallerThanParentCount += 1
					smallerThanParentCountPerSubclone += 1
		
			
			perSubcloneCounts[simulationDataIds[simulationDataSet]][vertexInd] = smallerThanParentCountPerSubclone
			
		if vertexCombinations == 0:
			smallerThanParentCounts.append(0)
			continue
		
		#The total score is out of how many combinations we see that the distance to the parent is larger than to other subclones. 
		smallerThanParentCounts.append(smallerThanParentCount / float(vertexCombinations))
		
		if (smallerThanParentCount / float(vertexCombinations)) > 0:
			
			print trees[simulationDataSet].edgeList
		
		
		
		#Then based on the alleles, compute the distance to the parent.
		
		#Then also compute the distance to all non-parents of this subclone
	
	return smallerThanParentCounts, perSubcloneCounts
		
#Also compute the ambiguity scores for each simulation round to determine if the subclone/parent distances correlate with failing to infer ambiguities correctly

def computeAmbiguousPositions():
	
	#We can generate the allele list with the event distances function
	kmin = 1 #the kmin and kmax used in the simulations
	kmax = 6 
	eventDistance = EventDistances(kmin, kmax)
	
	#get the allele list
	alleleList = eventDistance.alleleList
	
	#make sure that alleles are not duplicated
	LAFAndCombinations = dict()
	normalAlleles = Alleles(1,1)
	for allele in alleleList:
		AOccurrences = [m.start() for m in re.finditer('A', allele)]
		ACount = len(AOccurrences)
		BOccurrences = [m.start() for m in re.finditer('B', allele)]
		BCount = len(BOccurrences)
		
		alleleObj = Alleles(ACount, BCount)
		if BCount > ACount or BCount == ACount:
			alleleCombination = AlleleCombination([normalAlleles, alleleObj])
			
			for muIndex in range(0,101):
				
				LAF = alleleCombination.computeLAF(Mu(muIndex)) #compute the LAF that each combination would generate
				if LAF not in LAFAndCombinations.keys():
					#LAFAndCombinations[LAF] = []
					LAFAndCombinations[LAF] = 0
				
				#LAFAndCombinations[LAF].append((alleleObj.getAllelesAsString(), muIndex))
				LAFAndCombinations[LAF] += 1
			
	
	#print LAFAndCombinations
	#For every mu, we should check which LAF the combination with normal would generate
	
	#With this dictionary, we can check if a LAF has more than one solution. If true, then we can check and see if the position is correct or not. With that, we compute a score showing the
	#number of ambiguous positions that we were able to infer correctly. This score is a higher wow factor than 
	return LAFAndCombinations

def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines

#For each simulated dataset, we look at the real mu and the real A and compute the LAF. Does it correspond to a LAF that has multiple solutions?
#If yes, did we get the position correct?
#Perhaps, if no, what is the reason? Is it because we inferred another ambiguous solution, or did the method go wrong completely? (This can maybe be interesting as well)

def computeCorrectAmbiguityScore(LAFAndCombinations, simulationFolder):
	ambiguityScores = []
	ambiguities = []
	correctAmbiguityPositions = 0
	totalAmbiguousPositions = 0
	totalSize = 0
	#We need to read the actual A matrix values and also the mu
	normalAlleles = Alleles(1,1)
	#1. read the simulated A matrix
	
	allPerColAmbiguities = dict()
	
	for subdir, dirs, files in os.walk(simulationFolder):
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		for file in files:
			if re.match('RealA', file): #read the file and obtain the a matrix
				realAMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
			if re.match('RealMu', file): #also read the real mu
				realMu = collectErrorsFromFile(file, subdir)
				
			#Then load the inferred A and mu
			if re.match('EstimatedA', file): #read the file and obtain the a matrix
				estimatedAMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
			if re.match('EstimatedMu', file): #also read the real mu
				estimatedMu = collectErrorsFromFile(file, subdir)
			
		#Compute the LAF that each measurement in the real data would generate
		perColAmbiguityCount = dict()
		for row in range(0, realAMatrix.shape[0]):
			for col in range(0, realAMatrix.shape[1]):
				
				if col not in perColAmbiguityCount:
					perColAmbiguityCount[col] = realAMatrix.shape[0]
				
				totalSize += 1
				#generate allele object
				allele = realAMatrix[row][col]
				AOccurrences = [m.start() for m in re.finditer('A', allele)]
				ACount = len(AOccurrences)
				BOccurrences = [m.start() for m in re.finditer('B', allele)]
				BCount = len(BOccurrences)
				
				alleleObj = Alleles(ACount, BCount)
				alleleCombination = AlleleCombination([normalAlleles, alleleObj])
				
				#Compute the LAF this combination would generate
				muNormal = 1-(realMu[col])
				realMuObj = Mu(int(muNormal*100)) #this function only takes integer indices!
				realLAF = alleleCombination.computeLAF(realMuObj)
				
				#Check if this LAF is ambiguous y/n.
				ambiguousCount = LAFAndCombinations[realLAF]
				
				
				#If the ambiguous count > 1 and we are correct, we make a note of that.
				if ambiguousCount > 1:
					totalAmbiguousPositions += 1
					
					if realAMatrix[row][col] == estimatedAMatrix[row][col]:
						correctAmbiguityPositions += 1
						perColAmbiguityCount[col] -= 1 #Determine how many positions are wrong. 
				
		#Divide the ambiguity score by the total number of positions.
		#print correctAmbiguityPositions
		#print correctAmbiguityPositions / float(totalAmbiguousPositions)
		#print totalAmbiguousPositions / float(totalSize)
		#ambiguityScores.append(correctAmbiguityPositions / float(totalAmbiguousPositions)) #Reporting as % of ambiuguities
		#Reporting the ambiguity scores as the fraction of the total
		ambiguityScores.append(correctAmbiguityPositions / float(totalSize))
		
		ambiguities.append(totalAmbiguousPositions / float(totalSize))
		allPerColAmbiguities[subdir] = perColAmbiguityCount
	#Compute an average for every noise level.
	
	#convert to z-scores
	
	averageAmbiguityScore = sum(ambiguityScores) / float(len(ambiguityScores))
	averageAmbiguities = sum(ambiguities) / float(len(ambiguities))
	
	return [averageAmbiguities, averageAmbiguityScore, ambiguityScores, allPerColAmbiguities]

	


[trees, alleles, simulationDataIds] = readSimulationData(simulationFolder)
[distances, perSubcloneCounts] = computeDistances(trees, alleles, simulationDataIds)
LAFAndCombinations = computeAmbiguousPositions() #This is a pre-computed dictionary of all LAF and the associated number of combinations.
#Compute ambiguity scores
[averageAmbiguities, averageAmbiguityScore, ambiguityScores, allPerColAmbiguities] = computeCorrectAmbiguityScore(LAFAndCombinations, simulationFolder)

print allPerColAmbiguities


#Correlate
print np.corrcoef(distances,ambiguityScores)

#Do a different correlation, but then for each subclone individually.
#Make a vector format, with all simulation datasets and subclones thereof in the correct order. Then correlate the ambiguities and distance vectors.

perSubcloneCountsVector = []
correctAmbiguitiesVector = []
for simulationDataset in allPerColAmbiguities:
	
	subclonesAmbiguities = allPerColAmbiguities[simulationDataset]
	subclonesDistances = perSubcloneCounts[simulationDataset]
	
	for subclone in subclonesDistances: #use this format, because the ambiguities include the healthy subclone, for which we did not compute the distances. 
	
		perSubcloneCountsVector.append(subclonesDistances[subclone])
		correctAmbiguitiesVector.append(subclonesAmbiguities[subclone])

print len(perSubcloneCountsVector)
print len(correctAmbiguitiesVector)
	
print np.corrcoef(perSubcloneCountsVector, correctAmbiguitiesVector) #What we correlate is the number of times the distance to other subclones is larger than to the parent, and how many ambiguities are NOT resolved. 

# 
# import matplotlib.pyplot as plt
# 
# fig, ax = plt.subplots()
# 
# plt.boxplot(distances)
# plt.ylabel("Percentage of subclones with larger distance to other subclones")
# ax.set_xticklabels("")
# plt.savefig("DistanceToOtherSubclones.svg")



