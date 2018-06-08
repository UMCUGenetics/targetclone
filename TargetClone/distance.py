from alleles import Alleles
from editDistance import EditDistance
import numpy as np
from fst import FST
import settings

#Class to automatically compute the event distances between all possible alleles between kmin and kmax
#Pre-computation of all distances will save computational time, we just need to access the distances from a matrix. 
class EventDistances:
	
	distanceMatrix = None
	alleleList = None
	alleleListDict = None
	
	def __init__(self, kmin, kmax):
		self.generateAlleleList(kmin, kmax)
		self.generateEventDistanceMatrix(self.alleleList)
	
	#Make a list of all possible alleles between kmin and kmax
	def generateAlleleList(self, kmin, kmax):
		self.alleleList = []
		for k in range(kmin, kmax+1):
			#make Allele object
			for copyN in range(0,k+1):
				ACount = copyN
				BCount = k - copyN
				self.alleleList.append(Alleles(ACount, BCount).alleleString)
				
		#make a second, dictionary-based allele list for lookup speed
		self.alleleListDict = dict()
		addedVals = 0
		for k in range(kmin, kmax+1):
			#make Allele object
			for copyN in range(0,k+1):
				ACount = copyN
				BCount = k - copyN
				self.alleleListDict[Alleles(ACount, BCount).alleleString] = addedVals #remmber index
				addedVals += 1
	
	#Generate a matrix containing the distances between all pre-computed possible alleles	
	def generateEventDistanceMatrix(self, alleleList):
		self.distanceMatrix = np.empty([len(alleleList), len(alleleList)])
		editDistance = EditDistance()
		for allele in range(0, len(alleleList)):
			for otherAllele in range(0, len(alleleList)):
				#compute event distance
				dist = editDistance.editDistance(alleleList[allele], alleleList[otherAllele])
				self.distanceMatrix[allele, otherAllele] = dist
	
	#compute the normalized event distance between a parent and target
	def getEventDistance(self, parent, target):
		#get index from allele list
		#parentInd = self.alleleList.index(parent)
		#targetInd = self.alleleList.index(target)

		#try speed of dictionary-based lookup
		parentInd = self.alleleListDict[parent]
		targetInd = self.alleleListDict[target]
		
		return (1/float(self.distanceMatrix[parentInd, targetInd] + 1) )
	
	#Compute the event distance without normalization
	def getRawEventDistance(self, parent, target): #event distance without the immediate division. Duplicate code!
		#parentInd = self.alleleList.index(parent)
		#targetInd = self.alleleList.index(target)
		
		parentInd = self.alleleListDict[parent]
		targetInd = self.alleleListDict[target]
		
		return self.distanceMatrix[parentInd, targetInd]
	
	#given our data, which is a matrix of the most likely alleles per SNP position per subclone
	#compute a pairwise distance between subclones.
	#we use the assumption that if one allele has been lost, it can never be regained.
	#we are not certain about all positions due to noise. Thus, we can assign a penalty rather than a distance of infinite.
	#the data is a matrix with alleles, measurement positions in the rows and samples in the columns. The function computes a total pairwise distance (replaced by the FST)
	def computePairwiseAlleleDistance(self, data):
		penalty = settings.fst['lossOfMediumConfidentLohPenalty']
		distanceMatrix = np.zeros((data.shape[1], data.shape[1]))
		for sample1Ind in range(0, data.shape[1]):
			for sample2Ind in range(0, data.shape[1]):
				sample1Alleles = data[:,sample1Ind]
				sample2Alleles = data[:,sample2Ind]
				#compute the event distance between these alleles.
				totalDistance = 0
				
				for alleleInd in range(0, len(sample1Alleles)):
					if sample1Alleles[alleleInd] is None or sample2Alleles[alleleInd] is None:
						continue
					if type(sample1Alleles[alleleInd]) is not float and type(sample2Alleles[alleleInd]) is not float:
						eventDistance = self.getRawEventDistance(sample1Alleles[alleleInd].alleleString, sample2Alleles[alleleInd].alleleString)
						totalDistance += eventDistance
						#Check if alleles are lost that cannot be regained. If true, assign a penalty.
						#We can assign a higher penalty when more alleles need to be lost. 
						if sample1Alleles[alleleInd].ACount == 0 and sample2Alleles[alleleInd].ACount > 0:
							
							totalDistance += (penalty * sample2Alleles[alleleInd].ACount) #we multiply the penalty with how many As we have more than in the other sample
						if sample1Alleles[alleleInd].BCount == 0 and sample2Alleles[alleleInd].BCount > 0:
							totalDistance += (penalty * sample2Alleles[alleleInd].BCount)
						
				distanceMatrix[sample1Ind, sample2Ind] = totalDistance
		return distanceMatrix

#Handler to compute various kinds of pairwise distances		
class PairwiseDistance:
	
	#Use the FST to compute distances between copy numbers. Data is here a matrix with the copy numbers of measurements (rows) and samples (columns)
	def fstDistance(self, data):
		fst = FST()
		sampleCount = data.shape[1]
		dm = np.zeros((sampleCount, sampleCount))
		
		for i in range(0, sampleCount):
			for j in range(0, sampleCount):
				#extract the whole copy number profile, so the whole column
				sample1Profile = data[:,i]
				sample2Profile = data[:,j]

				dm[i,j] = fst.computeDistance(sample1Profile, sample2Profile)
		return dm
	
	#Compute the distances between subclones based on alleles using the FST. Data is a matrix with allele objects of measurements (rows) and samples (columns) 
	def fstAlleleDistance(self, data, samples): #this could have been more efficient when we pass objects
		fst = FST()
		sampleCount = data.shape[1]

		dm = np.zeros((sampleCount, sampleCount))
		messages = dict()
		for i in range(0, sampleCount):
			for j in range(0, sampleCount):
				#extract the whole copy number profile, so the whole column
				sample1Profile = data[:,i]
				sample2Profile = data[:,j]
				returnValue = fst.computeAlleleDistance(sample1Profile, sample2Profile, samples[i], samples[j])
				messages[(i,j)] = returnValue[0]

				dm[i,j] = returnValue[1]
		print "distances: "
		print dm
		return [messages, dm]
	
	#deprecated, quick way of computing absolute distance between all columns of a matrix
	def pdist(self, data):

		m = data.shape[1]
		n = data.shape[0]
		mask = np.isfinite(data)
		dm = np.zeros(m * (m - 1) // 2, dtype=float)
		dm2 = np.zeros((m,m))
		k = 0
		for i in range(0, m-1): #column 1
			for j in range(i + 1, m): #adjacent column
				#extact the whole column

				curr = np.logical_and(mask[:,i], mask[:,j])

				u = data[:,i][curr]
				v = data[:,j][curr]
				dm[k] = sum(abs(u - v))
				dm2[i][j] = dm[k]
				dm2[j][i] = dm[k]
				
				k += 1
		return dm2
	
#Handler for computing the distances given somatic variants.
#For each sample, the somatic variants must be made binary.
#Then, we extract the inferred C of the neighboring LAF
#If for either of these a copy is lost, then we tolerate the loss of a somatic variant. 
class SomaticVariantDistance:
	
	#Compute the distance between samples based on somatic variants
	def computeSampleDistance(self, samples):
		somVarDistanceMatrix = np.ones((len(samples),len(samples))) #Binary. Keep 1 on the diagonal as we will multiply this matrix with the distances based on c
		messages = dict()
		svMatrix = np.ones((len(samples[0].somaticVariants),len(samples)))
		for i in range(0, len(samples)): #skip sample 0, this is the healthy sample and can be the parent of everything by default. 
			#Extract the somatic variants and their indices
			#For each somatic variant, get the index. The index will be - subtraction for each somatic variant we process. 
			svMatrix[:,i] = samples[i].somaticVariants 
			for j in range(0, len(samples)):
				
				returnList = self.computeDistanceBetweenSomaticVariants(samples[i], samples[j],i,j)
				somVarDistanceMatrix[i][j] = returnList[1]
				distanceMessage = returnList[0]
				messages[(i,j)] = distanceMessage
				
		return [messages, somVarDistanceMatrix]
	
	#Compute the distance between samples based on somatic variants
	#The distance will be 0 if a relation is possible. If the relation is not possible (loss not tolerated), a penalty is assigned. 			
	def computeDistanceBetweenSomaticVariants(self, sample1, sample2,i,j):
		distance = 1
		somVarIndices = sample1.somaticVariantsInd
		messages = [] #This is the message that we will show in the tree for this relation between samples. 
		sample1Subtraction = 0
		for somVarInd in range(0, len(sample1.somaticVariants)):

			#If the somatic variant is NA, we don't know anything about this position, so skip it.
			if np.isnan(sample1.somaticVariants[somVarInd]) or np.isnan(sample2.somaticVariants[somVarInd]):
				continue
			
			#Extract the right measurements around this somatic variant (immediately adjacent, may not be enough evidence)
			noOfSomVarBefore = somVarInd+1
			#print "som var ind: ", somVarInd
			prevInd = (somVarIndices[somVarInd]-noOfSomVarBefore)-1
			nextInd = prevInd+1
			if prevInd+1 > len(sample1.measurements.measurements)-1:
				nextInd = prevInd
			if prevInd == -1:
				prevInd = 0
				
			# print sample1
			# print sample1.name
			# print len(sample1.bestCMu)
			# print prevInd
			# print nextInd

			cMuPrev1 = sample1.bestCMu[prevInd] 
			cMuNext1 = sample1.bestCMu[nextInd] #current ind
			
			#If both of these are None, we must not check it any further, as we do not have enough information to be confident.
			cMuPrev2 = sample2.bestCMu[prevInd]
			cMuNext2 = sample2.bestCMu[nextInd]
		
			#Either of these can be true, now BOTH need to be true for it to work!
			#If the previous is none, but the second one is not, we still have enough information
			copyLost = False
			none1 = False
			none2 = False
			if cMuPrev1 is not None and cMuPrev2 is not None:
				if cMuPrev2.c.c[1] - cMuPrev1.c.c[1] < 0: #we have lost a copy in the previous position
					copyLost = True
			else:
				none1 = True
			if cMuNext1 is not None and cMuNext2 is not None:
				if cMuNext2.c.c[1] - cMuNext1.c.c[1] < 0: #we have lost a copy in the next position
					copyLost = True
			else:
				none2 = True
			if none1 is True and none2 is True:
				copyLost = True #we do not have enough information to say if the variant is lost together with the allele or not, so we tolerate potential loss.
			
			
			#We have lost a copy and a somatic variant, append the message that the variant is lost. 
			if copyLost == True:
				if (sample1.somaticVariants[somVarInd] > 0 and sample2.somaticVariants[somVarInd] == 0):
					messages.append('- SNV ' + str(somVarInd + 1))
			if copyLost == False: #If we have not lost and copy and we do lose the somatic variant, then this is an impossible relation. 
				if (sample1.somaticVariants[somVarInd] > 0 and sample2.somaticVariants[somVarInd] == 0):
					messages.append('- SNV ' + str(somVarInd + 1))
					penalty = settings.fst['lossOfSNVPenalty'] #high penalty, loss could also be measurement noise, so do not use inf. 
					return [messages, penalty]
			if (sample1.somaticVariants[somVarInd] == 0 and sample2.somaticVariants[somVarInd] > 0):
				messages.append('+ SNV ' + str(somVarInd + 1))
			sample1Subtraction += 1
		
		return [messages, distance]
				
	
	
				
				
				