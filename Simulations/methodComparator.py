#This class will do the following:
#1. Read the output files from one simulation folder (noise 0.02)
#1.5 Prepare input for LiCHeE and MEDICC
#2. Run LiCHeE and MEDICC on the data from the simulation folder
import pickle
import os
import re
import sys
import numpy as np
import subprocess
from copy import deepcopy
import os.path

np.set_printoptions(threshold=np.nan)

#bad path usage
sys.path.insert(0, '/Users/mnieboer/Documents/Backups/Code/hpcTargetClone_3011/Simulations/')
sys.path.insert(0, '/Users/mnieboer/Documents/Backups/Code/hpcTargetClone_3011/TargetClone/')

from simulations_mixedSamples import SimulationProbabilities
from distance import PairwiseDistance
from alleles import Alleles
from segmentations import Segmentation

class MethodComparator:
	#from the given simulation data folder, go through each folder and obtain:
	#- The real copy numbers per sample (v)
	#- The real tree (v)
	#- The real somatic variants (v) + positions and chromosomes (v)
	
	def readSimulationData(self, simulationFolder):
		simulationClasses = dict()
		for subdir, dirs, files in os.walk(simulationFolder):
			if subdir == simulationFolder: #we are not interested in the root folder
				continue
			for file in files:
				if re.match('simulationData', file): #read the file and obtain the error
					
					pkl_file = open(subdir + "/" + file, 'rb')
					simulationData = pickle.load(pkl_file)
					pkl_file.close()
					simulationClasses[subdir] = simulationData #store per folder to obtain the right dir for writing later
					break
		return simulationClasses
	
	def getSNVCoordinates(self):
		simulationProbabilities = SimulationProbabilities()
		#- store the probabilities in an object (simulationprobabilities)
		simulationProbabilities.readFullProbabilityFile('/Users/mnieboer/Documents/Backups/Code/hpcTargetClone_3011/Simulations/simulationProbabilities.txt') #very quick and dirty
		#extract chromosomes and positions]
		chromosomes = []
		positions = []
		for variant in simulationProbabilities.somaticVariants:
			
			#remove the letter from the chromosomes
			
			chromosomeNum = variant.chromosome[:-1]

			chromosomes.append(chromosomeNum)
			positions.append(variant.position)

		return [chromosomes, positions]
			
	#To prepare our input for LiCHeE, we need to make the following format:
	#Chromosome	Position	Description	Normal	I1	I2	II1	II2	IV1	IV2	IV3	V1
	#Chromosome and position are not in the simulation data output, but likely only in the pkl
	#The sample names are also not there in the folder, only in the trees

	#The chromosomes and positions are for the AF, not for the SNV. These can be found in another file (simulationProbabilities)
	def prepareLiCHeEInput(self, simulationClasses, outfileName):
		
		#get the chromosomes and positions for the SNVs
		[chromosomes, positions] = self.getSNVCoordinates()
			
		#for each folder, get the required data
		for subdir in simulationClasses:
			print subdir
			f = open(subdir + '/' + outfileName, 'w')
			header = "Chromosome\tPosition\tDescription\tReference"
			#Add the sample names
			for sample in simulationClasses[subdir].samples:
				#if sample.name == 'Healthy':
				#	continue
				sampleName = sample.name
				header += "\t"
				header += sampleName
			f.write(header)
			f.write("\n")
			
			
			snvMatrix = simulationClasses[subdir].snvMatrix
			snvMatrix = np.loadtxt(subdir + '/SomVar_0.txt', dtype=float)
			
			
			
			#Format the input file and write it to the directory
			for positionInd in range(0, len(positions)):
				
				chromosome = chromosomes[positionInd]
				position = positions[positionInd]
				snvs = snvMatrix[positionInd,:]
				
				#reduce SNVs to 0.999 if 1 and 0.001 if 0
				newSNVs = deepcopy(snvs)
				for snvInd in range(0, len(snvs)):
					if snvs[snvInd] == 1:
						newSNVs[snvInd] = 0.99999
					else:
						newSNVs[snvInd] = 0.00001
					
				
				line = chromosome + "\t" + position + "\tDescription\t" + "0\t0\t" + "\t".join([str(snv) for snv in snvs]) #add the 0 for the healthy sample
				#line = chromosome + "\t" + position + "\tDescription\t" + "\t".join([str(snv) for snv in snvs]) #add the 0 for the healthy sample
				f.write(line)
				f.write("\n")
				
			f.close()
			
	def parseLiCHeEOutput(self, simulationClasses):
		
		distanceMap = dict()
		for subdir in simulationClasses:
			
			#Format inFile
			inFile = 'licheeInput.txt.dot'
			
			#Check if file exists, otherwise skip
			
			
			if not os.path.exists(subdir + '/' + inFile): #if there is no output file, skip it. 
				continue 
			#determine which are the samples
			sampleOrder = []
			sampleNames = []
			with open(subdir + '/' + inFile, "r") as f:
				
				lineNum = 0
				for line in f:
					
					line = line.strip('\r\n')
					#print line
					#check if the line contains a label numbers

					if re.search("labelloc", line) is not None:

						#split the sample names from these lines
						splitLine = re.split("\s", line)
						
						#Regex just doesn't seem to work for obtaining the sample name
						sampleNameLocation = splitLine[3] #the actual name of the sample, not the integer code
						splitSampleName = re.split('"', sampleNameLocation)
						realSampleName = splitSampleName[1]
						sampleOrder.append(realSampleName)
						
						sampleName = splitLine[0]
						
					
						if sampleName not in sampleNames:
							sampleNames.append(int(sampleName))

			
			edgeList = []
			sampleEdges = []
			with open(subdir + '/' + inFile, "r") as f:
				
				lineNum = 0
				for line in f:
					
					line = line.strip('\r\n')
					#print line
					#check if the line contains a label numbers
					if re.search("\-\>", line):
						
						#split the sample names from these lines
						sampleName = int(re.search("\-\>\s+(\d+)\[", line).group(1))
						
						parentName = int(re.search("(\d+)\s+\-\>\s+\d+\[", line).group(1))
						
						if sampleName in sampleNames:
							distance = float(re.search("label\=\"(\d+\.\d+)\"", line).group(1))
						else:
							distance = 0
						
						edge = (distance, parentName, sampleName)
						if edge not in edgeList:
							edgeList.append(edge)
						if sampleName in sampleNames:
							sampleEdges.append(edge)
	
			
			#The fw paths contains all of our paths
			#For each sample, we want the shortest path
			edgeConnector = EdgeConnector()
			distances = edgeConnector.start(edgeList, sampleEdges, sampleNames)
			distanceMap[subdir] = [distances, sampleOrder]
			
		return distanceMap
		
		
	def formatMEDICCInput(self, simulationClasses):
		
		#MEDICC needs a fasta-like output where the copy numbers are separated per sample per chromosome.
		#Each chromosome has a separate file, for the major and minor copy numbers.
		#We should read the original alleles file, and obtain the major and minor copy numbers fom there.
		#In the files, each sample is on a separate line, where there is also a line for the diploid.
		#Make sure to skip positions where the copy number is > 4.
		
		descFile = 'desc.txt' #in here store the chromosome file names
		
		#1. Make the chromosome file names
		chrFilesMajor = []
		chrFilesMinor = []
		chrs = range(1,23) #we leave out X and Y in this simulation
		majorSuffix = '_major.fasta'
		minorSuffix = '_minor.fasta'
		for chrom in chrs:
			
			chrFileMajor = 'chr' + str(chrom) + majorSuffix
			chrFileMinor = 'chr' + str(chrom) + minorSuffix
			chrFilesMajor.append(chrFileMajor)
			chrFilesMinor.append(chrFileMinor)
		
		
		
		
		#Go through the simulation data
		for subdir in simulationClasses:
			
			splitDir = re.split("\/", subdir)
			shortDir = 'MEDICCInput/' + splitDir[len(splitDir)-1]
			if not os.path.exists(shortDir):
				os.makedirs(shortDir)
				
			descFile = open(shortDir + '/desc.txt', 'w')	
			#make the desc file in each folder
			for chrom in range(0, len(chrs)):
				
				descFile.write('chr' + str(chrs[chrom]) + '\t' + chrFilesMajor[chrom] + '\t' + chrFilesMinor[chrom])
				descFile.write("\n")
				
			descFile.close()	
			
			fileInd = 0
			#Obtain the chromosome division to separate the files
			chromosomes = simulationClasses[subdir].chromosomes
			
			sampleNames = []
			
			for sample in simulationClasses[subdir].samples:
				sampleNames.append(sample.name)
			
			aMatrix = simulationClasses[subdir].aMatrix
			previousChromosome = chromosomes[0]
			currentCopyNumbers = [] #this should be a subset of the numpy matrix, but size cannot be dynamic. list of lists?
			new = True
			expectedArmCount = 2 #we expcet 2 arms, also create 2 data points in the file
			for chromosomeInd in range(0, len(chromosomes)):
				chromosome = chromosomes[chromosomeInd]
				
				
				#For each chromosome, obtain the sample copy numbers.
				#If the chromosome changes, write the copy numbers to the respective file
				#Then continue to the next file. 

				if chromosome != previousChromosome:
					
					#load a new file
					
					fMajor = open(shortDir + "/" + chrFilesMajor[fileInd], 'w')
					fMinor = open(shortDir + "/" + chrFilesMinor[fileInd], 'w')
					#write chromosomes to fasta file, each sample is a new line
					for sample in range(0, aMatrix.shape[1]):
						
						
						lengthAfterFilter = 0 #keep track of the positions that we skip because > 4
						
						if sampleNames[sample] != 'Healthy':
							
							header = ">" + sampleNames[sample]
							#Count the number of A and B in the alleles
							
							aCount = []
							bCount = []
							previousACount = -1
							previousBCount = -1
							for copyN in currentCopyNumbers[sample]:
		
								#Check if the value is different from the previous. If yes, save the current value. 
								if copyN.ACount != previousACount or copyN.BCount != previousBCount:
									
									if copyN.ACount > 4 or copyN.BCount > 4:
										continue
									
									if copyN.ACount > copyN.BCount: #switch major and minor when ACount is larger than BCount. 
										bCount.append(copyN.ACount)
										aCount.append(copyN.BCount)
									else:	
										aCount.append(copyN.ACount)
										bCount.append(copyN.BCount)
									
									lengthAfterFilter += 1 #we've added a position, keep track of that
										
									previousACount = copyN.ACount
									previousBCount = copyN.BCount
									
								else:
									continue
					
							
							
							
							chrStringMajor = "".join([str(x) for x in bCount])
							chrStringMinor = "".join([str(x) for x in aCount])
							
							if len(chrStringMajor) < expectedArmCount: #if we lack one posiiton, it means that both arms do not have a different value. Thus we can extend by itself. 
								chrStringMajor += chrStringMajor
							if len(chrStringMinor) < expectedArmCount: #if we lack one posiiton, it means that both arms do not have a different value. Thus we can extend by itself. 
								chrStringMinor += chrStringMinor

						else:
							header = ">diploid"
							chrStringMajor = "".join(['1']*expectedArmCount)
							chrStringMinor = chrStringMajor
							

						
						fMajor.write(header)
						fMajor.write("\n")
						fMajor.write(chrStringMajor)
						fMajor.write("\n")
						
						fMinor.write(header)
						fMinor.write("\n")
						fMinor.write(chrStringMinor)
						fMinor.write("\n")
						
						
						
					#clear out chromosomes
					fileInd += 1
					currentCopyNumbers = []
					fMajor.close()
					fMinor.close()
					previousChromosome = chromosome
					new = True
					
				else: #store the current copy numbers. We want this to be a list of lists (or a numpy matrix?) to be a list across the samples. 
					for sample in range(0, aMatrix.shape[1]):
						if new == True: #new indicates that we are looking at a new chromosome, and thus the copy numbers need to be filled. 
							currentCopyNumbers.append([])
						else:
							#Obtain the actual copy numbers (split major and minor)
							currentCopyNumbers[sample].append(aMatrix[chromosomeInd][sample])
			
					new = False
			
			#Repeat for the last chromosome
			fMajor = open(shortDir + "/" + chrFilesMajor[fileInd], 'w')
			fMinor = open(shortDir + "/" + chrFilesMinor[fileInd], 'w')
			#write chromosomes to fasta file, each sample is a new line
			for sample in range(0, aMatrix.shape[1]):
				
				
				lengthAfterFilter = 0 #keep track of the positions that we skip because > 4
				if sampleNames[sample] != 'Healthy':
					
					header = ">" + sampleNames[sample]
					#Count the number of A and B in the alleles
					
					aCount = []
					bCount = []
					
					#do a segmentation. Each chromosome has only 2 possible copy number values (1 per chromosome arm), these can be the same.
					#only save these max 2 values. This will greatly speed up MEDICC. 
					previousACount = -1
					previousBCount = -1
					for copyN in currentCopyNumbers[sample]:

						#Check if the value is different from the previous. If yes, save the current value. 
						if copyN.ACount != previousACount or copyN.BCount != previousBCount:
							
							if copyN.ACount > 4 or copyN.BCount > 4:
								continue
							
							if copyN.ACount > copyN.BCount: #switch major and minor when ACount is larger than BCount. 
								bCount.append(copyN.ACount)
								aCount.append(copyN.BCount)
							else:	
								aCount.append(copyN.ACount)
								bCount.append(copyN.BCount)
							
							lengthAfterFilter += 1 #we've added a position, keep track of that
								
							previousACount = copyN.ACount
							previousBCount = copyN.BCount
							
						else:
							continue
					
					
					chrStringMajor = "".join([str(x) for x in bCount])
					chrStringMinor = "".join([str(x) for x in aCount])
					
					if len(chrStringMajor) < expectedArmCount: #if we lack one posiiton, it means that both arms do not have a different value. Thus we can extend by itself. 
						chrStringMajor += chrStringMajor
					if len(chrStringMinor) < expectedArmCount: #if we lack one posiiton, it means that both arms do not have a different value. Thus we can extend by itself. 
						chrStringMinor += chrStringMinor

				else:
					header = ">diploid"
					chrStringMajor = "".join(['1']*expectedArmCount)
					chrStringMinor = chrStringMajor
					

				
				fMajor.write(header)
				fMajor.write("\n")
				fMajor.write(chrStringMajor)
				fMajor.write("\n")
				
				fMinor.write(header)
				fMinor.write("\n")
				fMinor.write(chrStringMinor)
				fMinor.write("\n")
			
		
	def parseMEDICCOutput(self, simulationClasses):
		
		outputSorted = dict() #for each simulation folder, store the distance matrix and also the order of the samples
		#1. obtain directory where medicc output is stored (will be the same as the output folder)
		for subdir in simulationClasses:
			
			splitDir = re.split("\/", subdir)
			shortDir = 'MEDICCInput/' + splitDir[len(splitDir)-1] + '/output/'
		
			#read the tree_final.dist for the distances
			inFile = shortDir + 'tree_final.dist'
			
			#Check if file exists, otherwise skip this one
			if not os.path.exists(inFile):
				print subdir
				continue
			
			
			sampleOrder = []
			sampleCount = 0
			diploidInd = 0
			with open(inFile, "r") as f:
				
				lineNum = 0
				for line in f:
					
					if lineNum == 0:
						
						sampleCount = int(line)
						lineNum += 1
						continue

					line = line.strip('\r\n')
					
					#find where the diploid sample is, skip this one later
					contents = re.split("\s+", line)
					
					#col 1 is the sample name, rest are the distances. We skip the last line, this is diploid (is it always the same?)
					
					sampleName = contents[0]
					
					if sampleName == 'diploid':
						diploidInd = lineNum
						sampleName = 'Healthy'
						
					sampleOrder.append(sampleName)

					lineNum += 1
					
			#then read the file again and use the distances to fillt he distance matrix
			distanceMatrix = np.zeros([sampleCount, sampleCount])
			
			
			with open(inFile, "r") as f:
				
				lineNum = 0
				sampleInd = 0
				for line in f:
					
					line = line.strip('\r\n')
					if lineNum == 0:	#skip the line with the sample count
						lineNum += 1
						continue
					lineNum += 1
					
					#obtain the distances from this line
					contents = re.split("\s+", line)
					#col 1 is the sample name, rest are the distances. We skip the last line, this is diploid (is it always the same?)
					
					#sampleName = contents[0]
					#if sampleName == 'diploid':
					#	continue
					
					matrixSampleInd = 0
					for contentInd in range(1, len(contents)):
						if contentInd == diploidInd:
							continue
						
						distance = float(contents[contentInd])
						
						distanceMatrix[sampleInd, matrixSampleInd] = distance
						matrixSampleInd += 1
						
					sampleInd += 1

						
			outputSorted[subdir] = [distanceMatrix, sampleOrder]
		
		return outputSorted
		
	#Get the distance matrix from the original trees	
	def parseGroundTruthDistanceMatrix(self, simulationClasses):
		segmentation = Segmentation()
		segmentation.setSegmentationFromFile('/Users/mnieboer/Documents/Backups/Code/hpcTargetClone_3011/TargetClone/InternalData/pq_segmentation.txt') #sorry
		#From the original simulation data, get the tree and generate a distance matrix
		distanceMatrices = dict()
		for subdir in simulationClasses:
			print subdir
			simulationClass = simulationClasses[subdir]
			
			aMatrix = simulationClass.aMatrix
			
			for sample in simulationClass.samples:
				sample.measurements.segmentation = segmentation
				#as the segmentation appears to be missing, we need to add this here
			
			
			distanceMatrix = PairwiseDistance().fstAlleleDistance(aMatrix, simulationClass.samples)
			print "subdir: ", subdir
			print distanceMatrix
			sampleOrder = []
			for sample in simulationClass.samples:
				sampleOrder.append(sample.name)
				
			distanceMatrices[subdir] = [distanceMatrix, sampleOrder]
		
		return distanceMatrices
	
	def parseTCOutputDistances(self, simulationClasses):
		#For this we need to read the actual output aMatrix files, from this we can again compute the distances based on all to all
		distanceMatrices = dict()
		for subdir in simulationClasses:

			simulationClass = simulationClasses[subdir]
			estimatedAStrMatrix = np.loadtxt(subdir + '/EstimatedA_0.txt', dtype='object')
			
			estimatedAMatrix = np.empty(estimatedAStrMatrix.shape, dtype='object')
			#convert the stirngs to allele objects
			for row in range(0,estimatedAStrMatrix.shape[0]):
				for col in range(0,estimatedAStrMatrix.shape[1]):
					#count the number of A's and B's in the string
					allele = estimatedAStrMatrix[row][col]
					AOccurrences = [m.start() for m in re.finditer('A', allele)]
					ACount = len(AOccurrences)
					BOccurrences = [m.start() for m in re.finditer('B', allele)]
					BCount = len(BOccurrences)
					
					alleleObj = Alleles(ACount, BCount)
					estimatedAMatrix[row][col] = alleleObj
					
			
			
			distanceMatrix = PairwiseDistance().fstAlleleDistance(estimatedAMatrix, simulationClass.samples)
			sampleOrder = []
			for sample in simulationClass.samples:
				sampleOrder.append(sample.name)
				
			distanceMatrices[subdir] = [distanceMatrix, sampleOrder]

		return distanceMatrices
		
		
	def rankSamplesByGroundTruth(self, simulationClasses, groundTruthDistances, tcEstimatedDistance, mediccDistances, liCHeEDistances):
		
		#loop through the ground truth, rank the samples by their distances.
		#Make two lists, one with the distances flattened, and one with the samples.
		
		rankedDistances = dict()
		#sort the lists by their distance, smallest distance first.
		treeErrors = []
		for subdir in groundTruthDistances:
			
			#Also check what the actual tree error is, does it correspond with the correlation that we get?
			treeError = np.loadtxt(subdir + '/treeError.txt', dtype='int')
			treeErrors.append(int(treeError))
			
			rankedDistances[subdir] = dict()
			distances = groundTruthDistances[subdir][0][1]
			sampleNames = groundTruthDistances[subdir][1]
			
			sampleCombinations = []
			distancesFlattened = []
			
			tcDistancesFlattened = []
			tcDistances = tcEstimatedDistance[subdir][0][1]
			#flatten sample names as: s1-s2, s1-s3, ...
			#Use a string such that we can easily find on which index the samples are
			#distances should be in the same format, keep as separate lists
			for sampleInd in range(0, len(sampleNames)):
				for sampleInd2 in range(sampleInd+1, len(sampleNames)):
					sampleCombination = sampleNames[sampleInd] + '-' + sampleNames[sampleInd2]
					sampleCombinations.append(sampleCombination)
					sampleCombination = sampleNames[sampleInd2] + '-' + sampleNames[sampleInd]
					sampleCombinations.append(sampleCombination)
					distancesFlattened.append(distances[sampleInd][sampleInd2])
					distancesFlattened.append(distances[sampleInd2][sampleInd])
					
					tcDistancesFlattened.append(tcDistances[sampleInd][sampleInd2])
					tcDistancesFlattened.append(tcDistances[sampleInd2][sampleInd])
					
			
			#sort both lists by the distance, ascending
			
			distancesFlattened = np.array(distancesFlattened)
			tcDistancesFlattened = np.array(tcDistancesFlattened)

			#get the indices of the sorted sample combinations
			
			sortedDistanceIndices = np.argsort(distancesFlattened)
			
			#then get the same ranked distances for the ground truth output
			rankedGTDistances = distancesFlattened[sortedDistanceIndices]
			
			#Sort the TC output in the same way (sample order is the same for the two)
			rankedTCDistances = tcDistancesFlattened[sortedDistanceIndices]


			#repeat for medicc output (here the sample order is different)
			#We need to check in which position the samples are and then store these at the right indices
			if subdir in mediccDistances:
				mediccSampleNames = mediccDistances[subdir][1]
				mediccEstimatedDistances = mediccDistances[subdir][0]
				rankedMediccDistances = deepcopy(rankedGTDistances)
				for sampleInd in range(0, len(mediccSampleNames)):
					for sampleInd2 in range(sampleInd+1, len(mediccSampleNames)):
						
						sampleCombination = mediccSampleNames[sampleInd] + '-' + mediccSampleNames[sampleInd2]
						#Find this sample combination in the sample combinations, check at which position it is in the array
				
						sampleIndex = sampleCombinations.index(sampleCombination)
						rankedMediccDistances[sampleIndex] = mediccEstimatedDistances[sampleInd][sampleInd2]
						
						#repeat for the reverse order
						sampleCombination = mediccSampleNames[sampleInd2] + '-' + mediccSampleNames[sampleInd]
						#Find this sample combination in the sample combinations, check at which position it is in the array
				
						sampleIndex = sampleCombinations.index(sampleCombination)
						rankedMediccDistances[sampleIndex] = mediccEstimatedDistances[sampleInd2][sampleInd]
				rankedDistances[subdir]['Medicc'] = rankedMediccDistances
		
			#Then finally lso repeat for the LiCHeE output. Some folders will be empty and have no output, but these can be skipped
			if subdir in liCHeEDistances:
				liCHeESampleNames = liCHeEDistances[subdir][1]
				liCHeEEstimatedDistances = liCHeEDistances[subdir][0]
				rankedLiCHeEDistances = deepcopy(rankedGTDistances)
				
				for sampleInd in range(0, len(liCHeESampleNames)):
					for sampleInd2 in range(sampleInd+1, len(liCHeESampleNames)):
						
						sampleCombination = liCHeESampleNames[sampleInd] + '-' + liCHeESampleNames[sampleInd2]
						#Find this sample combination in the sample combinations, check at which position it is in the array
	
						sampleIndex = sampleCombinations.index(sampleCombination)
						rankedLiCHeEDistances[sampleIndex] = liCHeEEstimatedDistances[sampleInd][sampleInd2]
						
						#repeat for the reverse order
						sampleCombination = liCHeESampleNames[sampleInd2] + '-' + liCHeESampleNames[sampleInd]
						#Find this sample combination in the sample combinations, check at which position it is in the array
	
						sampleIndex = sampleCombinations.index(sampleCombination)
						rankedLiCHeEDistances[sampleIndex] = liCHeEEstimatedDistances[sampleInd2][sampleInd]
				
				rankedDistances[subdir]['LiCHeE'] = rankedLiCHeEDistances
			
				
			rankedDistances[subdir]['GroundTruth'] = rankedGTDistances
			rankedDistances[subdir]['TargetClone'] = rankedTCDistances
		
		return rankedDistances
	
	def computeCorrelation(self, rankedDistances):
		tcCorrelationList = []
		mediccCorrelationList = []
		licheeCorrelationList = []
		correlation = dict()
		#compute correlation for each tool with the ground truth
		for subdir in rankedDistances:
			distances = rankedDistances[subdir]
			
			#sometimes the ground truth distancea and TC distances can be inf, we should appropriately mask these values or use a very large value.
			distances['GroundTruth'][np.where(distances['GroundTruth'] == float("inf"))] = 100000
			distances['TargetClone'][np.where(distances['TargetClone'] == float("inf"))] = 100000
			
			tcCorrelation = np.corrcoef(distances['GroundTruth'],distances['TargetClone'])
			
			correlation[subdir] = dict()
			correlation[subdir]['tcCorrelation'] = tcCorrelation[0][1]
			
			
			if 'Medicc' in distances:
				mediccCorrelation = np.corrcoef(distances['GroundTruth'],distances['Medicc'])
				correlation[subdir]['mediccCorrelation'] = mediccCorrelation[0][1]
				mediccCorrelationList.append(mediccCorrelation[0][1])
			
			if 'LiCHeE' in distances:
				licheeCorrelation = np.corrcoef(distances['GroundTruth'], distances['LiCHeE'])
				correlation[subdir]['licheeCorrelation'] = licheeCorrelation[0][1]
				licheeCorrelationList.append(licheeCorrelation[0][1])
				if licheeCorrelation[0][1] < 0:
					print str(subdir) + "LICHEE: " + str(licheeCorrelation[0][1])
					print str(subdir) + "MEDICC: " + str(mediccCorrelation[0][1])
					print str(subdir) + "TC: " + str(tcCorrelation[0][1])
					
			
			
			
			tcCorrelationList.append(tcCorrelation[0][1])
			
			
		
		# multiple box plots on one figure
		tcCorrelationList = np.array(tcCorrelationList)
		tcCorrelationList = tcCorrelationList[np.logical_not(np.isnan(tcCorrelationList))] #remove the nan values for testing
		mediccCorrelationList = np.array(mediccCorrelationList)
		licheeCorrelationList = np.array(licheeCorrelationList)
		print len(licheeCorrelationList)
		
		
		
		fig = figure()
		ax = axes()
		hold(True)
		boxplot(tcCorrelationList, positions = [1.5], widths=0.1)
		boxplot(mediccCorrelationList, positions = [2], widths=0.1)
		boxplot(licheeCorrelationList, positions = [2.5], widths=0.1)
		xticks([1.5, 2, 2.5], ['TargetClone', 'MEDICC', 'LiCHeE'], rotation=90)
		xlim(1, 3)
		ylim(-1.5,1.5)
		ylabel('Ranked correlation with ground truth')
		fig.tight_layout()
		#show()
		savefig('simulationMethodsCorrelation.svg')


	##ALTERNATIVE PLOTTING: for when we only want to compare the results for which LiCHeE was able to reconstruct a tree. 
	def computeCorrelationSubset(self, rankedDistances):
		tcCorrelationList = []
		mediccCorrelationList = []
		licheeCorrelationList = []
		correlation = dict()
		#compute correlation for each tool with the ground truth
		for subdir in rankedDistances:
			distances = rankedDistances[subdir]
			
			#sometimes the ground truth distancea and TC distances can be inf, we should appropriately mask these values or use a very large value.
			distances['GroundTruth'][np.where(distances['GroundTruth'] == float("inf"))] = 100000
			distances['TargetClone'][np.where(distances['TargetClone'] == float("inf"))] = 100000
			
			tcCorrelation = np.corrcoef(distances['GroundTruth'],distances['TargetClone'])
			
			correlation[subdir] = dict()

			
			if 'LiCHeE' in distances:
				licheeCorrelation = np.corrcoef(distances['GroundTruth'], distances['LiCHeE'])
				correlation[subdir]['licheeCorrelation'] = licheeCorrelation[0][1]
				licheeCorrelationList.append(licheeCorrelation[0][1])
				
				if 'Medicc' in distances: #only add MEDICC if there is also LiCHeE
					mediccCorrelation = np.corrcoef(distances['GroundTruth'],distances['Medicc'])
					correlation[subdir]['mediccCorrelation'] = mediccCorrelation[0][1]
					mediccCorrelationList.append(mediccCorrelation[0][1])
			
				#Only add TC when there is also LiCHeE
				correlation[subdir]['tcCorrelation'] = tcCorrelation[0][1]
				tcCorrelationList.append(tcCorrelation[0][1])
			
			
		
		# multiple box plots on one figure
		tcCorrelationList = np.array(tcCorrelationList)
		tcCorrelationList = tcCorrelationList[np.logical_not(np.isnan(tcCorrelationList))] #remove the nan values for testing
		mediccCorrelationList = np.array(mediccCorrelationList)
		licheeCorrelationList = np.array(licheeCorrelationList)
		
		
		
		fig = figure()
		ax = axes()
		hold(True)
		boxplot(tcCorrelationList, positions = [1.5], widths=0.1)
		boxplot(mediccCorrelationList, positions = [2], widths=0.1)
		boxplot(licheeCorrelationList, positions = [2.5], widths=0.1)
		xticks([1.5, 2, 2.5], ['TargetClone', 'MEDICC', 'LiCHeE'], rotation=90)
		xlim(1, 3)
		ylim(-1.5,1.5)
		ylabel('Ranked correlation with ground truth')
		fig.tight_layout()
		#show()
		savefig('simulationMethodsCorrelation_LiCHeESubset.svg')
		

class EdgeConnector:
	
	#Make a recursive function where we start at one edge and continue until we hit a leaf node. All these edges should be stored with their distances
	def findEdgeBW(self, edgeList, sampleNames, currentEdge, currentPath, distance, paths, pathDistances):
		
		for edge in edgeList:
			
			matchFound = False #keep track of if there is another edge, otherwise this path ends here and does not have other connections. 
			#Starting from one sample, we search bottom up for the next node that connects to it.
			#So the parent of the current sample needs to be found
			if edge[2] == currentEdge[1]:
				
				matchFound = True
				#extend our current path
				currentPath.append(edge[2])
				distance += edge[0]
				
				[paths, pathDistances] = self.findEdgeBW(edgeList, sampleNames, edge, currentPath, distance, paths, pathDistances)
				
				if currentPath not in paths:
					paths.append(currentPath)
					pathDistances.append(distance)
			
				#Actually, we need to check if there really is another path. If not, we also do not need to add it.
				#If the edge[1] is nowhere else in edge[2], then we add it.
				
				duplicate = False
				for edge2 in edgeList:
					if edge[1] == edge2[2]:
						duplicate = False
				
				if duplicate == False:
					
					duplicatePath = deepcopy(currentPath)
					duplicatePath.append(edge[1])
					if duplicatePath not in paths:
						paths.append(duplicatePath)
						pathDistances.append(distance) #the distance should be different?
						
			if matchFound == False:
				if [edge[2], edge[1]] not in paths:
					paths.append([edge[2], edge[1]])
					pathDistances.append(edge[0])
		return [paths, pathDistances]
					
	#If we have collected our backwards paths, we continue extending these forward until we end up at the samples again. 
	def findEdgeFW(self, edgeList, sampleNames, currentPath, distance, fwPaths, fwPathDistances):
		
		for edge in edgeList:
			
			#if we cannot find any more hits here, we need to continue searching for the next node the last node of the current path connects to
			if edge[1] == currentPath[(len(currentPath)-1)]:
				
				#again extend the path
				newPath = deepcopy(currentPath)
				newPath.append(edge[2])
				
				newDistance = distance + edge[0]
				[fwPaths, fwPathDistances, foundPath] = self.findEdgeFW(edgeList, sampleNames, newPath, newDistance, fwPaths, fwPathDistances)
				
				#print "current path now: ", currentPath
				if foundPath[(len(foundPath)-1)] in sampleNames and foundPath[0] != foundPath[(len(foundPath)-1)]:
					if foundPath not in fwPaths:
						fwPaths.append(foundPath) #only append the path when we are really done. 
						fwPathDistances.append(newDistance)
					
		return [fwPaths, fwPathDistances, currentPath]
	
	def start(self, edgeList, sampleEdges, sampleNames):
				
		paths = []
		fwPaths = []
		pathDistances = []
		fwPathDistances = []
		
		for sampleEdge in sampleEdges:
			
			currentPath = [sampleEdge[2]]
			[paths, pathDistances] = self.findEdgeBW(edgeList, sampleNames, sampleEdge, currentPath, sampleEdge[0], paths, pathDistances)
		
		for path in range(0, len(paths)):
			
			[fwPaths, fwPathDistances, path] = self.findEdgeFW(edgeList, sampleNames, paths[path], pathDistances[path], fwPaths, fwPathDistances)
		
		
		#Make unique sample names
		uniqueSampleNames = []
		for sampleName in sampleNames:
			if sampleName not in uniqueSampleNames:
				
				uniqueSampleNames.append(sampleName)
		
		#Make a matrix where the distance represents the shortest distance path for these two samples.
		import numpy as np
		distances = np.zeros([len(uniqueSampleNames), len(uniqueSampleNames)])
		for sampleName1Ind in range(0, len(uniqueSampleNames)):
			sampleName1 = uniqueSampleNames[sampleName1Ind]
			for sampleName2Ind in range(0, len(uniqueSampleNames)):
				sampleName2 = uniqueSampleNames[sampleName2Ind]
				currentDistance = float("inf")
				bestPathInd = 0
				for path in range(0, len(fwPaths)):
					if fwPaths[path][0] == sampleName1 and fwPaths[path][(len(fwPaths[path])-1)] == sampleName2:
						if fwPathDistances[path] < currentDistance:
							currentDistance = fwPathDistances[path]
							bestPathInd = path
				
				distances[sampleName1Ind][sampleName2Ind] = currentDistance		
				
		
		return distances		


				
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, xticks, ylabel	

methodComparator = MethodComparator()
#Obtain all simulations
print "reading input data..."
simulationClasses = methodComparator.readSimulationData('/Users/mnieboer/Documents/Backups/Code/hpcTargetClone_3011/Results/noise0RealVAF_mac/')
print "done"
#For each of these simulations, create a file that can be used as input to LiCHeE
#print "preparing input for LiCHeE"
#methodComparator.prepareLiCHeEInput(simulationClasses, 'licheeInput.txt')
#print "done"

#parse LiCHeE output to distance matrix
print "parsing LiCHeE output..."
liCHeEDistances = methodComparator.parseLiCHeEOutput(simulationClasses)
print "done"

#Format the input data for MEDICC
#print "formatting input for MEDICC..."
#methodComparator.formatMEDICCInput(simulationClasses)
#print "done"

#Parse the MEDICC output and generate a distance matrix
print "parsing medicc output..."
mediccDistances = methodComparator.parseMEDICCOutput(simulationClasses)
print "done"
# 
# #Read the original simulation data, obtain the distance matrix (truth)
print "parsing ground truth distances..."
groundTruthDistances = methodComparator.parseGroundTruthDistanceMatrix(simulationClasses)
print "done"
# 
# #Read the TargetClone output distances
print "parsing TC distances..."
tcEstimatedDistance = methodComparator.parseTCOutputDistances(simulationClasses)
print "done"

# 
# #Rank the samples by smallest to largest distance in the original simulation, keep sample order the same between tools
print "ranking distances..."
rankedDistances = methodComparator.rankSamplesByGroundTruth(simulationClasses, groundTruthDistances, tcEstimatedDistance, mediccDistances, liCHeEDistances) #very static and hardcoded for now, could be nicer but I won't bother
print "done"

# 
# #Compute a correlation, use this as a score of how good all tools are.
print "computing correlation..."
methodComparator.computeCorrelation(rankedDistances)
print "done"

#alternative for only plotting the results for which LiCHeE found a tree
# print "computing correlation for subset..."
# methodComparator.computeCorrelationSubset(rankedDistances)
# print "done"