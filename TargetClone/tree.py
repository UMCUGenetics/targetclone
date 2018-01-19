from collections import defaultdict, namedtuple
import numpy as np
from copy import deepcopy
from distance import PairwiseDistance
from fst import FST
import math

from zss import simple_distance, Node

#Describes all functions to make on a graph/tree object. 
class Graph:
	
	vertices = None
	edges = None
	graph = None
	edgeList = None #for easy searching
	edgeAnnotations = None #we also annotate edges with what happens. LOH, or loss of somatic variant
	#For each edge, we store the edge name and then the event that happened. 
	
	def __init__(self, vertices, edges, edgeList):
		self.setVertices(vertices)
		self.setEdges(edges)
		self.setGraph()
		self.edgeList = edgeList
		self.edgeAnnotations = dict()

	def setVertices(self, vertices):
		self.vertices = vertices
	
	#if the edges are changed, make sure that the edge list is updated as well 	
	def setEdges(self, edges):
		self.edges = edges
		if edges is None:
			self.edgeList = None
		else:
			self.edgeList = list(edges)
		self.graph = {'vertices' : self.vertices, 'edges' : self.edges}
		
	def setGraph(self, graph = None):
		if graph is None:
			self.graph = {'vertices' : self.vertices, 'edges' : self.edges}
		else:
			self.graph = graph
			
	def getGraph(self):
		return self.graph
	
	def addVertices(self, vertices):
		self.vertices.append(vertices)
		self.graph = {'vertices' : self.vertices, 'edges' : self.edges}
		
	def addEdgeAnnotation(self, edge, annotation):
		self.edgeAnnotations[edge] = annotation
	
	#sometimes it is easy to set the edges and their weight immediately from a distance matrix
	def setEdgesFromDistances(self, distances):
		
		edgeList = list()
		for i in range(0, distances.shape[0]):
			for j in range(0, distances.shape[1]):
				edgeList.append((distances[i][j], i, j))
		
		self.edges = set(edgeList)
		self.setGraph(None)
		self.edgeList = edgeList

	#given the name of a sample, search through the edges and find the parent of this sample
	def getParentByVertexId(self, vertexId):
		#Loop through the edges, third entry.
		for edge in self.edgeList:
			if edge[2] == vertexId:
				return edge[1]
		return 0 #the healthy sample/precursor is the root of the tree so that each subclone always has a parent.
	
	#remove specified edge from the tree
	def removeEdge(self, edge):
		edgeListCopy = deepcopy(self.edgeList)
		
		for selfEdge in range(0, len(self.edgeList)):
			if self.edgeList[selfEdge][1] == edge[1] and self.edgeList[selfEdge][2] == edge[2]:
				del edgeListCopy[selfEdge]
		
		self.setEdges(edgeListCopy) #and update the full graph after we edited its edges
		
	def addEdge(self, edge):
		self.edgeList.append(edge)
		#update the graph too
		self.graph = {'vertices' : self.vertices, 'edges' : set(self.edgeList)}
	
	#Function will check if the given graph is the same as this graph. The edges must be the same, and the weights too.  
	def compareIfGraphIsTheSame(self, graph):
		if self.getTotalWeight() == graph.getTotalWeight():
			#Check if the edges are the same, regardless of the weight. Sometimes some weights only differ slightly, this is not worth it. 
			selfEdges = []
			alternativeEdges = []
			for edge in range(0, len(self.edgeList)):
				selfEdges.append((self.edgeList[edge][1], self.edgeList[edge][2]))
				alternativeEdges.append((graph.edgeList[edge][1], graph.edgeList[edge][2]))
			if set(selfEdges) == set(alternativeEdges):
				return True
		return False
	
	#A tree is valid when we sum the introduction of a new somatic variant across all branches is at most 1 (infinite sites assumption).
	#Bad edges are the ones that are involved in getting a sum larger than 1.
	#These edges should be returned by the method.
	def checkIfTreeIsValid(self, variants): #returns true if the tree is valid, otherwise gets the next iteration. 
		#Check if this means if any somatic variant is lost
		#This goes wrong because some precursors are not added, but their identifier is increased. 
		adjacencyMatrix = np.zeros((len(self.graph['vertices']),len(self.graph['vertices'])))
		
		for edge in range(0, len(self.graph['edges'])):
			node = self.edgeList[edge][1]
			child = self.edgeList[edge][2]
			adjacencyMatrix[node][child] = 1

		badEdges = []
		
		presenceArray = np.ma.zeros((len(self.edgeList), variants.shape[0])) #for each branch we store if it introduces a new somatic variant compared to the parent or not. 
		
		#For every branch, get the somatic variants and the variants of the parent. Take the absolute distance and store this in the presence array.
		
		#If a variant is NA, we cannot say anything about it. 
		
		binaryVariants = np.zeros(variants.shape)
		binaryVariants[np.where(variants > 0)] = 1
		binaryVariants[np.where(np.isnan(variants))] = np.nan
		
		maskedVariants = np.ma.array(binaryVariants, mask = False)
		maskedVariants.mask[np.isnan(binaryVariants)] = True
		
		#look at the variants. If a variant is lost we do not count it as being harmful to the ISA as we do not know about re-gains from this.
		#If a variant is gained, we add a distance of 1 to the matrix. 
		for edgeInd in range(0,len(self.edgeList)):
			parent = self.edgeList[edgeInd][1]
			node = self.edgeList[edgeInd][2]
			
			#We need to keep in mind here that when a somatic variant is lost, this is not bad!!!
			#Also keep in mind that when we have no information on the parent, we can still say something about presence in the child. 
			distance = np.ma.array(np.zeros(len(maskedVariants[:,parent])), mask = False)
			
			for variant in range(0, len(maskedVariants[:,parent])): #we can probably do this in a faster way
				
				#If a somatic variant comes from NA, we really cannot say anything. 
				if maskedVariants.mask[variant][node] == True or maskedVariants.mask[variant][parent] == True: 
					distance[variant] = 0
					
				elif maskedVariants[variant][parent] > 0 and maskedVariants[variant][node] == 0:
					distance[variant] = 0 #we do not count this because a loss does not mean gained twice independently.
				else:
					distance[variant] = np.absolute(maskedVariants[variant][parent] - maskedVariants[variant][node])

			
			#store the distances in the presence array, this shows us exactly at which edges which variants are added (edges x variants (rxc))
			presenceArray[edgeInd,:] = distance

		#Sum across the columns to see which variants are introduced more than once
		summedPresence = np.sum(presenceArray, axis=0)
		
		#return the edges that are involved in having a summed presence larger than 1.
		#These are all edges that have a 1 in a column where the sum is larger than 1.
		badVarInd = np.where(summedPresence > 1)[0]
		badEdgeList = np.zeros(len(self.edgeList), dtype=int)
		for ind in badVarInd:
			involvedInd = np.where(presenceArray[:,ind] > 0)
			badEdgeList[involvedInd] = 1
			
		#obtain the edges in list format	
		badEdgeList = np.where(badEdgeList > 0)[0]
		
		#if there are no bad edges, we are done resolving the ISA. 
		if len(badEdgeList) < 1:
			return [None, None, None]
		#Make a new version of the presence array but then only with the involved edges.
		involvedPresence = presenceArray[badEdgeList,:]

		#check which edge needs to be removed based on which has the largest distance to the parent.
		worstEdgeInd = np.argwhere(involvedPresence.sum(axis=1) == np.max(involvedPresence.sum(axis=1)))
		
		worstEdgeList = [] #convert to a proper list
		for edge in worstEdgeInd:
			worstEdgeList.append(edge[0])
		
		#worstEdgeInd = np.argwhere(involvedPresence.sum(axis=1) == np.max(involvedPresence.sum(axis=1)))[0]
		#print "worst edges: ", worstEdgeInd
		#The edges with the highest distances are potential edges to remove. 
		#potentialEdgesToRemove = badEdgeList[worstEdgeInd]
		potentialEdgesToRemove = badEdgeList[worstEdgeList]
		
		#obtain the actual edges
		badEdges = [ self.edgeList[i] for i in potentialEdgesToRemove]
		allBadEdges = [ self.edgeList[i] for i in badEdgeList]
		
		#report the worst edges. 
		if badEdges is not None: 
			return [badEdges, allBadEdges, presenceArray]
		else:
			noneList = [None, None, None]
			return noneList
	
	#Report the total weight of the graph. 	
	def getTotalWeight(self):
		totalWeight = 0
		for edge in self.edgeList:
			totalWeight += edge[0]
		return totalWeight
	
	#Option to add one precursor node (unsampled subclone)
	def addPrecursorNode(self, originalGraph, originalSomaticVariants, samples, seenPrecursorSamples, dists):
		print "X"
		print dists
		#Make a backup of the original states so that we can reset
		fullGraph = deepcopy(originalGraph) #make a copy
		newSomaticVariants = deepcopy(originalSomaticVariants)
		
		#copy variants to a masked matrix so that we can ignore the NA positions
		binarySomVar = np.ma.array(newSomaticVariants, mask = False)
		binarySomVar[np.where(newSomaticVariants > 0)] = 1
		binarySomVar.mask[np.isnan(binarySomVar)] = True
		
		
		#we search for the sample that has the lowest number of somatic variants (column sum).
		#we first try to make a new subclone with the lowest number of variants and try to see if this one resolves the ISA. 
		colSum = np.sum(binarySomVar, axis=0)
		sortedSums = np.sort(colSum)
		sortedSumsInd = np.argsort(colSum)
	
		currentBest = 0
		
		for nextSumInd in range(0, len(sortedSums)):
			if sortedSums[nextSumInd] > currentBest and sortedSums[nextSumInd] not in seenPrecursorSamples:
				currentBest = sortedSumsInd[nextSumInd]
				seenPrecursorSamples.append(sortedSums[nextSumInd])
				break
		seenPrecursorSamples.append(currentBest)
		print "current best: ", currentBest
		if currentBest == 0:
			print "Could not resolve the ISA with precursors, reporting the best minimum distance tree instead"
			unresolved = True
			return [None, None, None]
		
		#we now know which somatic variants we need to copy.
		#if the column that we are interested in has a 1, we also set this to 1. If NA, we also set it to NA to be representative of the other clone. 
		allSomaticVariantsDup = np.zeros((newSomaticVariants.shape[0],newSomaticVariants.shape[1]+1))
		allSomaticVariantsDup[:,:-1] = newSomaticVariants

		for newVarInd in range(0, newSomaticVariants.shape[0]):
			if binarySomVar[newVarInd][currentBest] > 0:
				allSomaticVariantsDup[newVarInd][(newSomaticVariants.shape[1])] = 1
			elif binarySomVar[newVarInd][currentBest] == 0:
				allSomaticVariantsDup[newVarInd][(newSomaticVariants.shape[1])] = 0
			else:
				allSomaticVariantsDup[newVarInd][(newSomaticVariants.shape[1])] = np.nan
		
		#update the somatic variants matrix to include the unsampled subclone
		newSomaticVariants = allSomaticVariantsDup
		
		#Add new branches, all relations to the new subclone should be possible
		currentEdges = fullGraph.edgeList
		
		#Compute a distance based on somatic variants
		#precursors should also be allowed to be parents of each other!
		#if a somatic variant is lost, we consider the relation to not be possible (we have no information on the alleles of the newly introduced subclone), this is pure speculation, update if necessary
		#otherwise, simply copy the distance from the parent of the subclone that the variants were copied from. we assume that these subclones are relatively similar. 
		
		edgeMessages = dict()
		for sample in range(0, len(samples)):

			distance = newSomaticVariants[:,currentBest] - newSomaticVariants[:,sample]
			#if any somatic variant is lost, do not tolerate it. 

			dist = 0
			distR = 0
			#a precursor represents that the somatic variant was still there, thus it should not be lost. 
			if len(np.where(distance > 0)[0]) > 0: #if the current best loses the somatic variant compared to the sample, the distance is 1.
				
				dist = float("inf")
				#tolerate a loss if this is also tolerated by the previous parent. 
				
			if len(np.where(distance < 0)[0]) > 0:
				distR = float("inf")
			
			if sample == currentBest: #avoid mask
				dist = 0
				distR = 0
			#else:
				#dist = dists[currentBest][sample]
				#distR = dists[sample][currentBest]
				#dist = 0
				#distR = 0
			
			#only add an edge if the distance is not infinite (restricted)
			if dist != float("inf"):
				currentEdges.append((dist, len(samples), sample))
			if distR != float("inf"):
				currentEdges.append((distR, sample, len(samples)))
			
			#we are also interested in which somatic variants are lost/gained in the precursor. Add these annotations too!
			messages = []
			reverseMessages = []
			for variant in range(0, len(newSomaticVariants[:,currentBest])):
				if newSomaticVariants[variant][currentBest] == 0 and newSomaticVariants[variant][sample] > 0:
					messages.append('+ SNV' + str(variant+1))
					reverseMessages.append('- SNV' + str(variant+1)) #the messages work the other way around for sample to precursor
				if newSomaticVariants[variant][currentBest] > 0 and newSomaticVariants[variant][sample] == 0:
					messages.append('- SNV' + str(variant+1))
					reverseMessages.append('+ SNV' + str(variant+1))
			#are we adding this to the correct branch?
			print len(samples)
			print sample
			print reverseMessages
			fullGraph.edgeAnnotations[(len(samples), sample)] = messages
			fullGraph.edgeAnnotations[(sample, len(samples))] = reverseMessages
		
		#update the full graoh with the new information. 
		fullGraph.setEdges(currentEdges)
		fullGraph.addVertices(len(samples))
		
		
		#re-make the tree, now with the new edges and the additional somatic variants representing a precursor sample. 
		newEdges = SpanningArborescence().computeMinimumSpanningArborescence(fullGraph, newSomaticVariants)
		
		#back to the original algorithm, we then try to resolve the ISA with this new precursor introduced. We will get back here if the precursor did not help, then we select the next best subclone to copy and retry.
		return [newEdges, newSomaticVariants, seenPrecursorSamples, fullGraph]
	
	#The Zhang-Shasha algorithm requires a format with Nodes to compute an edit distance. So, we convert our tree to that format. 
	def convertTreeToEditDistanceFormat(self):
		nodeObjects = dict()
		
		#remove all weights and then make a set. This should help in having the same ordered list. If the order of the edges is not the same, then zss will mess up the edit distance somehow. 
		unweightedEdges = []
		for edge in self.edgeList:
			newEdge = (0, edge[1], edge[2])
			unweightedEdges.append(newEdge)
		
		unweightedEdges.sort(key=lambda tup: tup[1])
		
		for edge in set(unweightedEdges):
			#obtain the parent and child
			parent = edge[1]
			child = edge[2]
			
			#Create an object for the parent.
			#The only object that we will be appending to is the parent one
			if parent not in nodeObjects.keys():
				parentNode = Node(str(parent), [])
				nodeObjects[parent] = parentNode
				
				
			if child not in nodeObjects.keys():
				childNode = Node(str(child), [])
				nodeObjects[child] = childNode
			
			nodeObjects[parent].addkid(nodeObjects[child])
		
		return nodeObjects
	
	#compute the zhang-shasha edit distance for the current tree compared to the given tree. 
	def computeTreeEditDistance(self, otherTree):
		
		ownTreeEditDistFormat = self.convertTreeToEditDistanceFormat()
		otherTreeEditDistFormat = otherTree.convertTreeToEditDistanceFormat()
		
		#now compute the edit distance between these two
		if 'Healthy' in ownTreeEditDistFormat.keys():
			distance = simple_distance(ownTreeEditDistFormat['Healthy'], otherTreeEditDistFormat['Healthy']) #from the root we can reach all the nodes.
		else:
			distance = simple_distance(ownTreeEditDistFormat['Precursor'], otherTreeEditDistFormat['Precursor']) #from the root we can reach all the nodes.

		return distance

#Edmond's algorithm to infer a minimum spanning arborescence for a graph with possible edges. 	
class SpanningArborescence:
	
	def computeMinimumSpanningArborescence(self, graph, variants):
		#given this distance matrix, make a spanning arborescence
		arcs = self.defineArcs(graph)
		arborescence = self.min_spanning_arborescence(arcs, 0) #root the tree at sample 0, this can be changed depending on what we know (?)
		if arborescence is False:
			return False
		#Report the arborescence in the above graph format (just stick with it even though it's another conversion.... what if we need to change back, then it's bad to
		#change our interface format)
		edges = []
		for arc in arborescence:
			edges.append((arborescence[arc][1], arborescence[arc][0], arborescence[arc][2]))
		return set(edges)
		
	
	#Given a distance matrix, define the Arc tuples for the subclones
	def defineArcs(self, graph):
		#use the input graph to define the arcs. 
		
		Arc = namedtuple('Arc', ('head', 'weight', 'tail'))
		arcList = []
		for edge in graph.edges: 
			arcList.append(Arc(edge[1], edge[0], edge[2]))


		return arcList
	
	#Run the algorithm to find an mst
	def min_spanning_arborescence(self, arcs, sink):
		good_arcs = []
		quotient_map = {arc.tail: arc.tail for arc in arcs}
		quotient_map[sink] = sink
		while True:
			min_arc_by_tail_rep = {}
			successor_rep = {}
			for arc in arcs:
				if arc.tail == sink:
					continue
				tail_rep = quotient_map[arc.tail]
				head_rep = quotient_map[arc.head]
				if tail_rep == head_rep:
					continue
				if tail_rep not in min_arc_by_tail_rep or min_arc_by_tail_rep[tail_rep].weight >= arc.weight: #set node to tail if the weight is larger?? If we use equals, we get more linear structures
					min_arc_by_tail_rep[tail_rep] = arc
					successor_rep[tail_rep] = head_rep
			cycle_reps = self.find_cycle(successor_rep, sink)
			if cycle_reps is False:
				return False
			if cycle_reps is None:
				good_arcs.extend(min_arc_by_tail_rep.values())
				return self.spanning_arborescence(good_arcs, sink)
			good_arcs.extend(min_arc_by_tail_rep[cycle_rep] for cycle_rep in cycle_reps)
			cycle_rep_set = set(cycle_reps)
			cycle_rep = cycle_rep_set.pop()
			quotient_map = {node: cycle_rep if node_rep in cycle_rep_set else node_rep for node, node_rep in quotient_map.items()}
	
	
	def find_cycle(self, successor, sink):
		visited = {sink}
		for node in successor:
			cycle = []
			while node not in visited:
				visited.add(node)
				cycle.append(node)
				if node not in successor:
					return False
				node = successor[node] #check if the node is there. If it is not, then our tree is not possible and we need to return this. 
			if node in cycle:
				return cycle[cycle.index(node):]
		return None
	
	
	def spanning_arborescence(self, arcs, sink):
		arcs_by_head = defaultdict(list)
		for arc in arcs:
			if arc.tail == sink:
				continue
			arcs_by_head[arc.head].append(arc)
		solution_arc_by_tail = {}
		stack = arcs_by_head[sink]
		while stack:
			arc = stack.pop()
			if arc.tail in solution_arc_by_tail:
				continue
			solution_arc_by_tail[arc.tail] = arc
			stack.extend(arcs_by_head[arc.tail])
		return solution_arc_by_tail


		