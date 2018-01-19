import numpy as np

class SimulationErrorHandler:
	
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
	
	def computeMuErrorFromVectors(self, realMu, eMu):
		totalMuDifference = 0
		for muInd in range(0, len(realMu)):
			
			
			muDifference = np.absolute(realMu[muInd] - eMu[muInd])
			totalMuDifference += muDifference
		
		muError = totalMuDifference / len(realMu)
		return muError

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
		#print "the best tree: ", bestTree.edgeList
		
		#If the best tree does not have any edges, we do not need to compare.
		if len(bestTree.edgeList) < 1:
			print "tree does not have edges, entire tree is incorrect"
			#Get the number of nodes in the real tree
			numberOfEdges = len(realTree.edgeList)
			treeDistance = numberOfEdges
		else:
			#compare the best tree to the edges of our tree.
			treeDistance = bestTree.computeTreeEditDistance(realTree)

		return treeDistance #use this for now.
	
	
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
				
				#This should not be possible, fix it later!
				if totalACount < 0:
					totalACount = 0
				if totalBCount < 0:
					totalBCount = 0
				countSum = totalACount + totalBCount
				minCount = min([totalACount] + [totalBCount])
				
				if countSum == 0: #this should not happen!!!!
					realLaf = 0
				else:
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
	