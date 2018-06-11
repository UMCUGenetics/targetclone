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
			print subdir
			
			if subdir == simulationFolder: #we are not interested in the root folder
				continue
			for file in files:
				if re.match('cError', file): #read the file and obtain the error
					cErrors += collectErrorsFromFile(file, subdir)
				if re.match('aError', file): #read the file and obtain the error
					aErrors += collectErrorsFromFile(file, subdir)
				if re.match('muError', file): #read the file and obtain the error
					muErrors += collectErrorsFromFile(file, subdir)
				if re.match('treeError', file): #read the file and obtain the error
					treeErrors += collectErrorsFromFile(file, subdir)

					
					
				#Repeat for the permutation errors
				if re.match('pCError', file): #read the file and obtain the error
					pCErrors += collectErrorsFromFile(file, subdir)
				if re.match('pAError', file): #read the file and obtain the error
					pAErrors += collectErrorsFromFile(file, subdir)
				if re.match('pMuError', file): #read the file and obtain the error
					pMuErrors += collectErrorsFromFile(file, subdir)
				if re.match('pTreeError', file): #read the file and obtain the error
					pTreeErrors += collectErrorsFromFile(file, subdir)
			
			
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
		groupedEuclideanErrors[noiseLevel] = euclideanErrors
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors]

readDataIncludingPermutations(motherFolder, snpNums)

#3. Make a plot with the error on the y axis, and the number of SNPs on the x axis. (4 plots per data type that we infer)


#4. Another plot that we can make is show how much performance is gained with the horizontal dependency on the y axis, and the average distance between SNPs on the x axis (all in 1 plot is possible)



