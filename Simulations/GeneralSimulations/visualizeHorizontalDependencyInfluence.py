"""
	The purpose of this script is to take as input a set of simulation results. One where the horizontal dependency exists and one where the LAF measurements are shuffled within a sample.
	The number of SNPs increases, meaning that we show exactly what the influence is of the horizontal dependency as the number of measurements decreases or increases.
	
"""
import sys
import os
import re

#1. Define the folders containing the data that we want to visualize
snpNums = [100]
motherFolder = sys.argv[1]

#2. Read the errors from these folders for the normal case and corresponding shuffling


def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines

def readDataIncludingPermutations(dataFolder, noiseLevels):
	
	groupedCErrors = dict()
	groupedAErrors = dict()
	groupedMuErrors = dict()
	groupedTreeErrors = dict()
	groupedAmbiguityErrors = dict()
	groupedPCErrors = dict()
	groupedPAErrors = dict()
	groupedPMuErrors = dict()
	groupedPTreeErrors = dict()
	groupedPAmbiguityErrors = dict()
	
	groupedEuclideanErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + '/snps_' + str(noiseLevel)
		
		#Read all the errors into one list for this noise level
		cErrors = []
		aErrors = []
		muErrors = []
		treeErrors = []
		ambiguityErrors = []
		ambiguityCorrectedErrors = []
		pCErrors = []
		pAErrors = []
		pMuErrors = []
		pTreeErrors = []
		pAmbiguityErrors = []
		pAmbiguityCorrectedErrors = []
		
		for subdir, dirs, files in os.walk(simulationFolder):
			
			
			if subdir == simulationFolder: #we are not interested in the root folder
				continue
			
			if re.search('horizontalShuffle', subdir) is not None:
				print subdir
				for file in files:
					if re.match('cError', file): #read the file and obtain the error
						pCErrors += collectErrorsFromFile(file, subdir)
						print "read: ", pCErrors
					if re.match('aError', file): #read the file and obtain the error
						pAErrors += collectErrorsFromFile(file, subdir)
					if re.match('muError', file): #read the file and obtain the error
						pMuErrors += collectErrorsFromFile(file, subdir)
					if re.match('treeError', file): #read the file and obtain the error
						pTreeErrors += collectErrorsFromFile(file, subdir)
			else:
				print "non shuffle: ", subdir
				for file in files:
					if re.match('cError', file): #read the file and obtain the error
						cErrors += collectErrorsFromFile(file, subdir)
					if re.match('aError', file): #read the file and obtain the error
						aErrors += collectErrorsFromFile(file, subdir)
					if re.match('muError', file): #read the file and obtain the error
						muErrors += collectErrorsFromFile(file, subdir)
					if re.match('treeError', file): #read the file and obtain the error
						treeErrors += collectErrorsFromFile(file, subdir)
						
			
		#Gather the data per noise level
		groupedCErrors[noiseLevel] = cErrors
		groupedAErrors[noiseLevel] = aErrors
		groupedMuErrors[noiseLevel] = muErrors
		groupedTreeErrors[noiseLevel] = treeErrors
		groupedAmbiguityErrors[noiseLevel] = ambiguityCorrectedErrors
		groupedPCErrors[noiseLevel] = pCErrors
		groupedPAErrors[noiseLevel] = pAErrors
		groupedPMuErrors[noiseLevel] = pMuErrors
		groupedPTreeErrors[noiseLevel] = pTreeErrors
		groupedPAmbiguityErrors[noiseLevel] = pAmbiguityCorrectedErrors
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors]

def averageData(dictionary, type):
	
	#Average the values per key in the dictionary
	#Make sure that the order matches
	averagedData = []
	noiseLevels = []
	for key in dictionary.keys():
		
		if len(dictionary[key]) < 1:
			return [] #sometimes we do not read certain data, skip it 
		
		noiseLevels.append(key)
		values = dictionary[key]
		newValues = values
		if type is not 'T':
			
			newValues = []
			for value in values:
				if value < 1: #somehow this happens for the permutations? Does it represent a float version of NA? We can ignore it, we have 10 permutations per. 
					newValues.append(value)
		
		average = sum(newValues)/float(len(newValues))
		averagedData.append(average)
		
	
	#Sort the averagedData
	sortedNoiseLevels, sortedAveragedData = (list(t) for t in zip(*sorted(zip(noiseLevels, averagedData))))
	
	return sortedAveragedData


#Get the raw errors for the horizontal shuffle and the normal case
[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors] = readDataIncludingPermutations(motherFolder, snpNums)

#Compute an average of the errors
print "C"
print averageData(groupedCErrors, 'C')
print averageData(groupedPCErrors, 'C')
print "A"
print averageData(groupedAErrors, 'C')
print averageData(groupedPAErrors, 'C')
print "Mu"
print averageData(groupedMuErrors, 'C')
print averageData(groupedPMuErrors, 'C')
print "T"
print averageData(groupedTreeErrors, 'T')
print averageData(groupedPTreeErrors, 'T')

#Compute the standard deviation of the error (add later)


#3. Make a plot with the error on the y axis, and the number of SNPs on the x axis. (4 plots per data type that we infer)


#4. Another plot that we can make is show how much performance is gained with the horizontal dependency on the y axis, and the average distance between SNPs on the x axis (all in 1 plot is possible)



