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

simulationFolder = "Results/generic_noise0.02/"

#From the simulation data, get the alleles and the tree from a file.
def readSimulationData(dataFolder):
	
	trees = []
	alleles = []
	
	for subdir, dirs, files in os.walk(simulationFolder):
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
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
			
		
	return trees, alleles

#For every subclone, compute the distance to its parent and to the other subclones. 
def computeDistances(trees, alleles):
	
	smallerThanParentCounts = []
	
	for simulationDataSet in range(0, len(trees)):
		
		print simulationDataSet
		
		currentTree = trees[simulationDataSet]
		currentAlleles = alleles[simulationDataSet]
		
		
		#For every subclone in the vertices list, determine who the parent is.
		for vertexInd in range(0, len(currentTree.vertices)):
			vertex = currentTree.vertices[vertexInd]
			parent = ""
			for edge in currentTree.edgeList:
				if edge[2] == vertex:
					parent = edge[1]
			
			if parent == "": #skip the subclones that have no parent
				continue
			
			# if parent == "Healthy": #it is difficult to use the healthy subclone in this analysis, because of the large number of events introduced in round 1 in particular. 
			# 	continue
					
			parentInd = currentTree.vertices.index(parent)
			
			#Compute the allelic event distance from the subclone to the parent
			subcloneAlleles = currentAlleles[:,vertexInd]
			parentAlleles = currentAlleles[:,parentInd]
			
			#Compute the event distance
			distanceToParent = 0
			for position in range(0, len(subcloneAlleles)):
				eventDistance = EventDistances(0, 6).getEventDistance(parentAlleles[position], subcloneAlleles[position]) #should use the settings
				distanceToParent += eventDistance

			#Now do the same for the non-parent subclones, and compute how often the distance is smaller than the distance to the parent.
			
			smallerThanParentCount = 0
			otherSubcloneCount = 0
			for vertexInd2 in range(vertexInd, len(currentTree.vertices)):
				vertex2 = currentTree.vertices[vertexInd2]
				if vertex2 == parent:
					continue
				
				otherSubcloneCount += 1
				
				totalDistance = 0
				otherSubcloneAlleles = currentAlleles[:,vertexInd2]
				for position in range(0, len(subcloneAlleles)):
					eventDistance = EventDistances(0, 6).getEventDistance(otherSubcloneAlleles[position], subcloneAlleles[position]) #should use the settings
					totalDistance += eventDistance
				if totalDistance < distanceToParent:
					smallerThanParentCount += 1
			
			#Make this a percentage. 	
			smallerThanParentCounts.append(smallerThanParentCount / float(otherSubcloneCount))
			
		
		#Then based on the alleles, compute the distance to the parent.
		
		#Then also compute the distance to all non-parents of this subclone
	
	print smallerThanParentCounts	
	return smallerThanParentCounts
		
		
	


[trees, alleles] = readSimulationData(simulationFolder)
computeDistances(trees, alleles)
