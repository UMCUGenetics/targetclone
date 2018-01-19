import numpy as np
from laf import LAF

import settings #file with all tool settings

#This should formally be an interface, with LAFSample as a specific type of sample with specific measurements.
#More formally, each sample can also consist of one or more subclone objects. These subclone objects will have A/C/Mu etc.
class Sample:
	#all the characteristics of a sample
	name = None
	measurements = None
	afMeasurements = None
	somaticVariants = None 
	somaticVariantsInd = None #we could make another class here for somatic variants, but that causes too much overhead.
	C = None
	A = None
	Mu = None
	sampleFile = None
	parent = None
	bestCMu = None
	originalBestCMu = None
	originalA = None #we need this to store the parent-of-parent information. It won't work otherwise as we will update. 
	
	def __init__(self, sampleFile, parent):
		self.sampleFile = sampleFile
		self.parent = parent
		self.bestCMu = []
		self.originalBestCMu = []
		self.A = []
	
	#The sample names are provided in the data file as column header
	def setName(self, name):
		self.name = name
	
	def setParent(self, parent): #We can use this to update the parent of a sample after each iteration. 
		self.parent = parent
		if self.parent is None:
			self.parent = self #If we have a healthy sample, this does not formally have a parent. We assume that the parent of the healthy sample is itself (also healthy), because we have no other information. 
		
	#The measurements sometimes contain NA strings, convert these to nan to be consistent in nan types
	def correctMeasurements(self, content):	
		content = [x.strip('\r\n') for x in content]
		
		#Remember the NA positions, add these back later or let the model handle it correctly!. This is important, because otherwise we are not comparing the right laf positions later!
		for c in range(0, len(content)): 
			if content[c] == 'NA':
				content[c] = np.nan
			
		#convert to numpy, this is faster for processing numerical arrays 
		measurements = np.asarray(content, dtype='float')
		return measurements
	
	#Obtain separate lists with SNP measurements and SNV measurements (indicated by 'N' and 'Y')
	def splitSomaticVariantsAndAfMeasurements(self, measurements, somaticVariantsIndication):
		afMeasurements = []
		somaticVariants = []
		for measurement in range(0, len(measurements)):
			if somaticVariantsIndication[measurement] == 'N':
				afMeasurements.append(measurements[measurement])
			else:
				
				somaticVariants.append(measurements[measurement])
				
		return [afMeasurements, somaticVariants]
	
	#Obtains columns with information from the data input file and set these as properties in the sample class			
	def setDataFromRawMeasurements(self, sampleName, measurements, chromosomes, starts, somaticVariantsIndication):
		#Obtain the header and remove it.
		measurements = list(measurements)
		correctedMeasurements = self.correctMeasurements(measurements)
		
		#Also remove the header for the other columns
		chromosomes = list(chromosomes)
		starts = list(starts)
		somaticVariantsIndication = list(somaticVariantsIndication)
		
		#Here filter the somatic variants for the ones that should be excluded (only here the positions still match)
		[filteredMeasurements, filteredChromosomes, filteredStarts, filteredSomaticVariantsIndication] = self.filterExcludedSNVs(correctedMeasurements, chromosomes, starts, somaticVariantsIndication)
		
		#Extract the LAF and somatic variants
		splitMeasurements = self.splitSomaticVariantsAndAfMeasurements(filteredMeasurements, filteredSomaticVariantsIndication)
		laf = LAF(splitMeasurements[0], filteredChromosomes, filteredStarts, filteredStarts) #convert to the right type of measurements
		self.measurements = laf
		self.afMeasurements = splitMeasurements[0]
		self.setName(sampleName)
		
		#Check here if the somatic variant has a low AF. If true, do not include it in further analyses!
		afLowerBound = float(settings.general['snvLowerBound'])
		filteredSvMeasurements = []
		for svMeasurement in splitMeasurements[1]:
			if np.isnan(svMeasurement) == True:
				filteredSvMeasurements.append(np.nan)
			elif svMeasurement > afLowerBound:
				filteredSvMeasurements.append(svMeasurement)
			else:
				filteredSvMeasurements.append(0) #could also be nan, but we count it as not detected.
		
		self.somaticVariants = filteredSvMeasurements
		#Obtain the indices of the somatic variants among the measurements, we use this to check if adjacent positions show LOH or not
		self.somaticVariantsInd = [i for i,v in enumerate(filteredSomaticVariantsIndication) if v == 'Y']
		
		#also remove the filtered out somatic variants entries from the chromosome/start information
		#do this in reverse to ensure that the indices do not switch
		for somVarIndex in sorted(self.somaticVariantsInd, reverse=True):
			del filteredChromosomes[somVarIndex]
			del filteredStarts[somVarIndex]
	
	#Provide the start positions of the somatic variants and their chromosomes
	#Read the exclusion file and see which positions match
	def filterExcludedSNVs(self, measurements, chromosomes, starts, snvIndications):
		
		filteredMeasurements = []
		filteredChromosomes = []
		filteredStarts = []
		filteredSNVIndications = []
		
		badPositions = []
		
		#Read the file with excluded positions
		exclusionFile = settings.files['excludedSNVs']
		
		with open(exclusionFile, "r") as inFile:
			lineNum = 0
			for line in inFile:
				if lineNum == 0: #skip the header
					lineNum += 1
					continue
				
				#check the data in this line.
				#Split the line by tab
				line = line.strip('\r\n')
				
				splitLine = line.split('\t')
				chromosome = splitLine[0]
				start = splitLine[1]
				badPositions.append((chromosome,start))
				
				lineNum += 1
		
		#loop over the positions. Does it match with any SNV that we wish to exclude?
		for position in range(0, len(measurements)):
			match = False
			for badPosition in range(0, len(badPositions)):
				#if the position is a somatic variant, the chromosome matches and the start matches, we have a match and exclude it. Otherwise, append
				
				if snvIndications[position] == 'Y' and chromosomes[position] == badPositions[badPosition][0] and starts[position] == badPositions[badPosition][1]:
					
					match = True
				
			if match == False:
				filteredMeasurements.append(measurements[position])
				filteredChromosomes.append(chromosomes[position])
				filteredStarts.append(starts[position])
				filteredSNVIndications.append(snvIndications[position])				
					
		
		return [filteredMeasurements, filteredChromosomes, filteredStarts, filteredSNVIndications]
	
	#Correct the measurements to set NA correctly, then add to the sample properties. 	
	def setSomaticVariantsFromRawMeasurements(self, measurements, indices):
		correctedSomVar = self.correctMeasurements(measurements)
		correctedIndices = self.correctMeasurements(indices)
		self.somaticVariants = correctedSomVar
		self.somaticVariantsInd = correctedIndices.astype(int)