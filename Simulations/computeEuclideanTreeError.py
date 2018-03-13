import os
import re
import numpy as np
import ast
from tree import Graph
from tree import SpanningArborescence

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

def computeTreeError(lafMatrix, snvMatrix, realTree):
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
