import sys
sys.path.insert(0, '../../TargetClone/')

import os
import re
import numpy as np
import ast
from tree import Graph
from tree import SpanningArborescence
from fst import FST
from sample import Sample
from laf import LAF
import simulationSettings
from segmentations import Segmentation
from alleles import Alleles
from distance import SomaticVariantDistance
from c import C
from simulation_generic import SimulationProbabilities
from dummyCMu import DummyCMu

from simulationErrors import SimulationErrorHandler
from scipy.spatial import distance

def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(line)

	text_file.close()
	return floatLines


#Provide a child to start the search from
#Then go up the tree, there should be only one path, until we reach the root of the tree.
def isParentChildRelation(edgeList, parent, child):
	
	print "looking for parent: ", parent
	print "looking for child: ", child 
	
	#First obtain the parent of the child in the tree to initiate the search
	
	#Keep searching starting from a new edge until we have found the parent that we are looking for or until we reach the root		
	parentFound = False
	root = False
	newChild = child
	while parentFound == False and root == False:
		
		if len(edgeList) < 1:
			root = True
		childPresent = False
		for edge in edgeList:
			
			if newChild == "Healthy":
				root = True #we can't continue after this
			
			if edge[1] == parent and edge[2] == newChild:
				#print "found parent: ", edge[1]
				parentFound = True


			#print "edge child: ", edge[2]
			#print "new child: ", newChild
			if edge[2] == newChild:
				childPresent = True
				newChild = edge[1] #Start from the parent of this child
				
				#print "setting new child to: ", newChild
		
		if childPresent == False: #Make sure to also quit when trees are truncated and the child that we are looking for is not there in the tree
			root = True
				
					
	
	print parentFound
	
	return parentFound
		
#Ancestry swap error
def computeAncestrySwapError(realTree, inferredTree):
	"""
		The purpose of this function is to check how often pairs of parents and children are swapped compared to the real ground truth tree. 
	"""
	
	print "edges: ", realTree.edgeList
	print "edges: ", inferredTree.edgeList
	
	allSamples = realTree.vertices
	
	samplePairs = []
	for sampleInd in range(0, len(allSamples)):
		for sampleInd2 in range(sampleInd+1, len(allSamples)):
			
			samplePairs.append([allSamples[sampleInd], allSamples[sampleInd2]])
			samplePairs.append([allSamples[sampleInd2], allSamples[sampleInd]]) #Also keep the reverse

	
	#Go through the edges in the real tree, and decide for each pair which sample is the parent and which is the child. If these are in different branches, we should not count these.
	
	
	parentChildRelationStatusRealTree = []
	parentChildRelationStatusInferredTree = []
	for samplePairInd in range(0, len(samplePairs)):
		
		samplePair = samplePairs[samplePairInd]
		print samplePair
		
		#Look through the edges of the real tree.
		#Search for the child of the pair.
		#Then continue up the branch recursively until the root is found. If we do not find the parent, skip it.
		#If we do find it as a parent, then we have found a parent-child relation.
		#Repeat this for each combination of sample pairs.
		#Afterwards, we repeat these steps for the inferred tree, and compare how often a parent-child relation is now swapped (we can see this by the order of the sample pairs)
		
		relationshipStatus = isParentChildRelation(realTree.edgeList, samplePair[0], samplePair[1]) #check if the first sample is the parent of the second in the tree. 
		parentChildRelationStatusRealTree.append(relationshipStatus)
		
		#Repeat for the inferred tree
		relationshipStatus = isParentChildRelation(inferredTree.edgeList, samplePair[0], samplePair[1]) #check if the first sample is the parent of the second in the tree. 
		parentChildRelationStatusInferredTree.append(relationshipStatus)
	
	#For every 2 indices, if the first is 1 for the real tree, but the second is 1 for the inferred tree, relationships are swapped.
	
	print "real tree statuses: ", parentChildRelationStatusRealTree
	print "inferred tree statuses: ", parentChildRelationStatusInferredTree
	
	ancestrySwapCountAbsentInInferred = 0
	ancestrySwapCountPresentInInferred = 0
	
	
	for samplePairInd in range(0, len(parentChildRelationStatusRealTree)):
		
		if parentChildRelationStatusRealTree[samplePairInd] == True and parentChildRelationStatusInferredTree[samplePairInd] == False:
			
			ancestrySwapCountAbsentInInferred += 1 #This is a relation that is present in the real tree, but not in the inferred tree.
		
		if parentChildRelationStatusRealTree[samplePairInd] == False and parentChildRelationStatusInferredTree[samplePairInd] == True:
			ancestrySwapCountPresentInInferred += 1 #This relation is not in the real tree, but newly added in the inferred tree. 
			
	
	print "absent in inferred error: ", ancestrySwapCountAbsentInInferred
	print "present in inferred error: ", ancestrySwapCountPresentInInferred
	
	print "done computing swap errors"	
	
	return [ancestrySwapCountAbsentInInferred, ancestrySwapCountPresentInInferred]


def computeEuclideanDistanceBetweenSamples(sample1Laf, sample2Laf, sample1Snvs, sample2Snvs):
	
	#compute the distance row-wise
	lafEuclidean = distance.euclidean(sample1Laf, sample2Laf)
	snvEuclidean = distance.euclidean(sample1Snvs, sample2Snvs)
	snvEuclidean = 0
	#lafEuclidean = 0

	return lafEuclidean + snvEuclidean

def generateInitialTree(dists, vertexNames):
	fullGraph = Graph(vertexNames, None, None)
	fullGraph.setEdgesFromDistances(dists)
	
	return fullGraph

def computeMST(fullGraph, vertexNames):
	newGraph = Graph(vertexNames, set(), [])
	newEdges = SpanningArborescence().computeMinimumSpanningArborescence(fullGraph, None)
	newVertices = vertexNames
	newGraph = Graph(newVertices, None, None)
	newGraph.setEdges(newEdges)
	
	#Update the node names in the edges
	updatedEdges = []
	for edge in newGraph.edgeList:
		newEdge = (edge[0], vertexNames[edge[1]], vertexNames[edge[2]])
		updatedEdges.append(newEdge)
	newGraph.edgeList = updatedEdges
	newGraph.edges = set(updatedEdges)
	return newGraph


def computeEuclideanTreeError(lafMatrix, snvMatrix, realTree):
	sampleNum = lafMatrix.shape[1]
	
	#Compute the distance pairwise between samples
	distanceMatrix = np.empty([sampleNum,sampleNum], dtype=float)
	
	for sample1 in range(0, sampleNum):
		for sample2 in range(0, sampleNum):
			euclideanDist = computeEuclideanDistanceBetweenSamples(lafMatrix[:,sample1], lafMatrix[:,sample2], snvMatrix[:,sample1], snvMatrix[:,sample2])
			
			distanceMatrix[sample1,sample2] = euclideanDist
	
	#Compute the MST
	fullGraph = generateInitialTree(distanceMatrix, realTree.vertices)
	mst = computeMST(fullGraph, realTree.vertices)
	simulationErrorHandler = SimulationErrorHandler()
	treeScore = simulationErrorHandler.computeTreeError([mst], realTree)
	return treeScore

def computeCTreeError(cMatrix, realTree):
	sampleNum = cMatrix.shape[1]
	
	#Compute the distance pairwise between samples
	distanceMatrix = np.empty([sampleNum,sampleNum], dtype=float)
	
	for sample1 in range(0, sampleNum):
		for sample2 in range(0, sampleNum):
			
			#The distance can be computed for the entire column at once using the FST
			dist = FST().computeDistance(cMatrix[:,sample1], cMatrix[:,sample2])
			
			distanceMatrix[sample1,sample2] = dist
	
	
	#Compute the MST
	fullGraph = generateInitialTree(distanceMatrix, realTree.vertices)
	mst = computeMST(fullGraph, realTree.vertices)
	simulationErrorHandler = SimulationErrorHandler()
	treeScore = simulationErrorHandler.computeTreeError([mst], realTree)
	return treeScore

#This is a bit annoying but we need the positions of the SNVs which have not been stored in a file
def obtainSomaticVariantIndices():
	simulationProbabilities = SimulationProbabilities()
	#- store the probabilities in an object (simulationprobabilities)
	simulationProbabilityFile = simulationSettings.files['simulationProbabilityFile']
	simulationProbabilities.readFullProbabilityFile(simulationProbabilityFile)

	[chromosomes, positions, segmentation, chromosomeArms] = parseReferenceFile()
	offset = 0
	variantIndices = []
	for variant in simulationProbabilities.somaticVariants:
		
		position = variant.position
		chromosome = variant.chromosome

		for measurementPosition in range(0, len(positions)-1):
			
			if str(chromosome) == str(chromosomeArms[measurementPosition]): #the chromosome needs to be the same
				
				#for all the variants within the SNP range		
				if int(position) > int(positions[measurementPosition]) and int(position) < int(positions[measurementPosition + 1]):
					variantIndex = measurementPosition + offset
					variantIndices.append(variantIndex)
					offset += 1
				#Situation where the somatic variant comes before the SNP measurements
				if str(chromosome) != str(chromosomeArms[measurementPosition-1]) and int(position) <= int(positions[measurementPosition]):
					variantIndex = measurementPosition + offset
					variantIndices.append(variantIndex)
					offset += 1
				#Situation where the somatic variant comes after the SNP measurements
				if str(chromosome) != str(chromosomeArms[measurementPosition+1]) and int(position) >= int(positions[measurementPosition]):
					variantIndex = measurementPosition + offset
					variantIndices.append(variantIndex)
					offset += 1
	
	return [chromosomes, positions, variantIndices]
	

def computeSNVTreeError(snvMatrix, cMatrix, lafMatrix, realTree, variantIndices, chromosomes, positions):
	sampleNum = snvMatrix.shape[1]
	
	cObjMatrix = np.empty(cMatrix.shape, dtype=object)
	for row in range(0, cMatrix.shape[0]):
		for col in range(0, cMatrix.shape[1]):
			currentC = cMatrix[row][col]
			cObj = C([2,int(currentC)], []) #empty vector to stop initialization of allele combinations
			dummyCMu = DummyCMu()
			dummyCMu.c = cObj
			cObjMatrix[row][col] = dummyCMu
	
	#Compute the distance pairwise between samples
	distanceMatrix = np.empty([sampleNum,sampleNum], dtype=float)
	
	for sample1 in range(0, sampleNum):
		for sample2 in range(0, sampleNum):
			
			#Make the sample objects. These now need somatic variants and a CMu
			sample1Obj = Sample(None, None)
			sample1Obj.bestCMu = cObjMatrix[:,sample1]
			#the dummy c mu is actually a list of dummy c mu's, so we need to make one for each c
			#dummyCMu = DummyCMu()
			#dummyCMu.c = cObjMatrix[:,sample1]
			#sample1Obj.bestCMu = dummyCMu
			sample1Obj.somaticVariants = snvMatrix[:,sample1]
			sample1Obj.somaticVariantsInd = variantIndices
			sample1Obj.measurements = LAF(lafMatrix[:,sample1], chromosomes, positions, positions)
			sample2Obj = Sample(None, None)
			#dummyCMu = DummyCMu()
			#dummyCMu.c = cObjMatrix[:,sample2]
			sample2Obj.bestCMu = cObjMatrix[:,sample2]
			sample2Obj.somaticVariants = snvMatrix[:,sample2]
			sample2Obj.somaticVariantsInd = variantIndices
			sample2Obj.measurements = LAF(lafMatrix[:,sample2], chromosomes, positions, positions)
			#The distance can be computed for the entire column at once using the FST
			[messages, dist] = SomaticVariantDistance().computeDistanceBetweenSomaticVariants(sample1Obj, sample2Obj, sample1, sample2)
			
			
			distanceMatrix[sample1,sample2] = dist
	
	#Compute the MST
	fullGraph = generateInitialTree(distanceMatrix, realTree.vertices)
	mst = computeMST(fullGraph, realTree.vertices)
	simulationErrorHandler = SimulationErrorHandler()
	treeScore = simulationErrorHandler.computeTreeError([mst], realTree)
	return treeScore

#we need the chromosome information in order for this to work. 
def parseReferenceFile():
	#Read the segmentation, use this to annotate chromosome positions with p/q arm.
	chromosomes = []
	positions = []
	chromosomeArms = []
	referenceFile = simulationSettings.files['referenceFile'] #In reference.txt we have stored the selected positions from T3209 that are not 0 in the reference sample. We have the chromosome numbers and positions. From the
	#measured AF values we estimate a SD that we use to generate a noise distribution from which we sample AF measurements. We generate new measurements at these specified positions. 
	segmentationFile = simulationSettings.files['segmentationFile']
	
	segmentation = Segmentation()
	segmentation.setSegmentationFromFile(segmentationFile)
	
	#Read the reference file
	with open(referenceFile, "r") as inFile:
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
			
			chromosomes.append(chromosome)
			positions.append(position)
			referenceMeasurements.append(float(splitLine[2]))
			chromosomeArmName = segmentation.getSegmentName(chromosome, position, position)
			chromosomeArms.append(chromosomeArmName)
			lineInd += 1
	
	return [chromosomes, positions, segmentation, chromosomeArms]
		
		
	

def computeATreeError(aMatrix, lafMatrix, afMatrix, realTree, chromosomes, positions):
	sampleNum = aMatrix.shape[1]
	
	aObjMatrix = np.empty(aMatrix.shape, dtype=object)
	#Convert the a matrix to an actual allele matrix
	for row in range(0, aMatrix.shape[0]):
		for col in range(0, aMatrix.shape[1]):
			allele = aMatrix[row][col]
			AOccurrences = [m.start() for m in re.finditer('A', allele)]
			ACount = len(AOccurrences)
			BOccurrences = [m.start() for m in re.finditer('B', allele)]
			BCount = len(BOccurrences)
			
			alleleObj = Alleles(ACount, BCount)
			aObjMatrix[row][col] = alleleObj
	
	#Compute the distance pairwise between samples
	distanceMatrix = np.empty([sampleNum,sampleNum], dtype=float)
	
	for sample1 in range(0, sampleNum):
		for sample2 in range(0, sampleNum):
			#make a dummy sample object for the FST function
			sample1Obj = Sample(None, None)
			sample1Obj.measurements = LAF(lafMatrix[:,sample1], chromosomes, positions, positions)
			sample1Obj.measurements.segmentation = segmentation
			sample1Obj.afMeasurements = afMatrix[:,sample1]
			sample2Obj = Sample(None, None)
			sample2Obj.measurements = LAF(lafMatrix[:,sample2], chromosomes, positions, positions)
			sample2Obj.measurements.segmentation = segmentation
			sample2Obj.afMeasurements = afMatrix[:,sample2]
			
			#The distance can be computed for the entire column at once using the FST
			[messages, dist] = FST().computeAlleleDistance(aObjMatrix[:,sample1], aObjMatrix[:,sample2], sample1Obj, sample2Obj)
			distanceMatrix[sample1,sample2] = dist
	#print distanceMatrix
	#exit()
	#Compute the MST
	fullGraph = generateInitialTree(distanceMatrix, realTree.vertices)
	mst = computeMST(fullGraph, realTree.vertices)
	simulationErrorHandler = SimulationErrorHandler()
	treeScore = simulationErrorHandler.computeTreeError([mst], realTree)
	return treeScore

