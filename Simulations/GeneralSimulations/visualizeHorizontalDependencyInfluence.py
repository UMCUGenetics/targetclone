"""
	The purpose of this script is to take as input a set of simulation results. One where the horizontal dependency exists and one where the LAF measurements are shuffled within a sample.
	The number of SNPs increases, meaning that we show exactly what the influence is of the horizontal dependency as the number of measurements decreases or increases.
	
"""
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import computeTreeErrorOtherMetrics

#1. Define the folders containing the data that we want to visualize
motherFolder = sys.argv[1]

snpNums = [100, 500, 1000, 5000, 10000, 50000] #use this for now because the folder names were changed to non-numbers
#snpNums = [500]
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
	groupedAncestrySwapErrors = dict()
	groupedAmbiguityErrors = dict()
	groupedPCErrors = dict()
	groupedPAErrors = dict()
	groupedPMuErrors = dict()
	groupedPTreeErrors = dict()
	groupedPAmbiguityErrors = dict()
	groupedPAncestrySwapErrors = dict()
	
	groupedEuclideanErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + '/snps_' + str(noiseLevel)
		
		#Read all the errors into one list for this noise level
		cErrors = []
		aErrors = []
		muErrors = []
		treeErrors = []
		ancestrySwapErrors = []
		ambiguityErrors = []
		ambiguityCorrectedErrors = []
		pCErrors = []
		pAErrors = []
		pMuErrors = []
		pTreeErrors = []
		pAmbiguityErrors = []
		pAmbiguityCorrectedErrors = []
		pAncestrySwapErrors = []
		
		for subdir, dirs, files in os.walk(simulationFolder):
			
			
			if subdir == simulationFolder: #we are not interested in the root folder
				continue
			
			if re.search('horizontalShuffle', subdir) is not None:
				for file in files:
					if re.match('cError', file): #read the file and obtain the error
						pCErrors += collectErrorsFromFile(file, subdir)
					if re.match('aError', file): #read the file and obtain the error
						pAErrors += collectErrorsFromFile(file, subdir)
					if re.match('muError', file): #read the file and obtain the error
						pMuErrors += collectErrorsFromFile(file, subdir)
					if re.match('treeError', file): #read the file and obtain the error
						pTreeErrors += collectErrorsFromFile(file, subdir)
					if re.match('RealTrees', file): #read the file and obtain the error
						stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
						tree = eval(stringDict)
						realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
						treeSizes.append(len(realTree.edgeList))
					
					if re.match('EstimatedTrees', file): #read the file and obtain the error
						stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
						tree = eval(stringDict)
						inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
						
				
				#Compute the ancestry swap error
				[ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)

				summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
				pAncestrySwapErrors.append(summedError / float(noOfSamplePairs))	
			
			else:
				for file in files:
					if re.match('cError', file): #read the file and obtain the error
						cErrors += collectErrorsFromFile(file, subdir)
					if re.match('aError', file): #read the file and obtain the error
						aErrors += collectErrorsFromFile(file, subdir)
					if re.match('muError', file): #read the file and obtain the error
						muErrors += collectErrorsFromFile(file, subdir)
					if re.match('treeError', file): #read the file and obtain the error
						treeErrors += collectErrorsFromFile(file, subdir)
					if re.match('RealTrees', file): #read the file and obtain the error
						stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
						tree = eval(stringDict)
						realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
						treeSizes.append(len(realTree.edgeList))
					
					if re.match('EstimatedTrees', file): #read the file and obtain the error
						stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
						tree = eval(stringDict)
						inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
						
				
				#Compute the ancestry swap error
				[ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)

				summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
				ancestrySwapErrors.append(summedError / float(noOfSamplePairs))	
			
		#Gather the data per noise level
		groupedCErrors[noiseLevel] = cErrors
		groupedAErrors[noiseLevel] = aErrors
		groupedMuErrors[noiseLevel] = muErrors
		groupedTreeErrors[noiseLevel] = treeErrors
		groupedAncestrySwapErrors[noiseLevel] = ancestrySwapErrors
		groupedAmbiguityErrors[noiseLevel] = ambiguityCorrectedErrors
		groupedPCErrors[noiseLevel] = pCErrors
		groupedPAErrors[noiseLevel] = pAErrors
		groupedPMuErrors[noiseLevel] = pMuErrors
		groupedPTreeErrors[noiseLevel] = pTreeErrors
		groupedPAncestrySwapErrors[noiseLevel] = pAncestrySwapErrors
		groupedPAmbiguityErrors[noiseLevel] = pAmbiguityCorrectedErrors
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestrySwapErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, groupedPAncestrySwapErrors]

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

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t.ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def obtainStandardDeviations(groupedErrors, averagedError):
	
	#Compute standard deviations, see if it is better
	groupedAboveStd = []
	groupedBelowStd = []
	
	#It is more interesting to show the quantiles rather than the standard deviation
	for noiseLevelInd in range(0, len(groupedErrors.keys())):
		noiseValues = groupedErrors[groupedErrors.keys()[noiseLevelInd]]
		[mean, lower, upper] = mean_confidence_interval(noiseValues)
		q1 = np.std(noiseValues)
		q3 = np.std(noiseValues)
		groupedAboveStd.append(upper)
		groupedBelowStd.append(lower)
		
	sortedKeys, sortedBelow = zip(*sorted(zip(groupedErrors.keys(), groupedBelowStd)))
	sortedKeys, sortedAbove = zip(*sorted(zip(groupedErrors.keys(), groupedAboveStd)))
	return [sortedAbove, sortedBelow]

#Get the raw errors for the horizontal shuffle and the normal case
[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestrySwapErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, groupedPAncestrySwapErrors] = readDataIncludingPermutations(motherFolder, snpNums)

#Compute an average of the errors
print "C"
averagedCErrors = averageData(groupedCErrors, 'C')
averagedPCErrors = averageData(groupedPCErrors, 'C')
print "A"
averagedAErrors = averageData(groupedAErrors, 'C')
averagedPAErrors = averageData(groupedPAErrors, 'C')
print "Mu"
averagedMuErrors = averageData(groupedMuErrors, 'C')
averagedPMuErrors = averageData(groupedPMuErrors, 'C')
print "T"
averagedTreeErrors = averageData(groupedTreeErrors, 'T')
averagedPTreeErrors = averageData(groupedPTreeErrors, 'T')
print averagedTreeErrors
print averagedPTreeErrors

averagedAncestrySwapErrors = averageData(groupedAncestrySwapErrors, 'C')
averagedPAncestrySwapErrors = averageData(groupedPAncestrySwapErrors, 'C')

#Compute the standard deviation of the error (add later)
[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)
[groupedAboveStdAncestry, groupedBelowStdAncestry] = obtainStandardDeviations(groupedAncestrySwapErrors, averagedAncestrySwapErrors)

#Standard deviations for the permutation run
[groupedAboveStdCP, groupedBelowStdCP] = obtainStandardDeviations(groupedPCErrors, averagedPCErrors)
[groupedAboveStdAP, groupedBelowStdAP] = obtainStandardDeviations(groupedPAErrors, averagedPAErrors)
[groupedAboveStdMuP, groupedBelowStdMuP] = obtainStandardDeviations(groupedPMuErrors, averagedPMuErrors)
[groupedAboveStdTP, groupedBelowStdTP] = obtainStandardDeviations(groupedPTreeErrors, averagedPTreeErrors)
[groupedAboveStdAncestryP, groupedBelowStdAncestryP] = obtainStandardDeviations(groupedPAncestrySwapErrors, averagedPAncestrySwapErrors)


#3. Make a plot with the error on the y axis, and the number of SNPs on the x axis. (4 plots per data type that we infer)

def plotHorizontalDependencyInfluence(errors, pErrors, aboveStd, belowStd, aboveStdP, belowStdP, snpNums, plotType, title, colInd):
	
	#Take the average and standard deviations as input
	
	#Make a plot with error on the y axis
	#Number of SNPs on the x axis.
	#For each run (with and without) we have 2 points
	colors = ['#1f77b4', '#ff7f0e', '#d62728', '#9467bd', '#2ca02c', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	plt.figure()
	ax = plt.gca()
	
	legendLines = []
	
	correctedBelowStd = []
	for std in range(0, len(belowStd)):
		newStd = belowStd[std]
		if (errors[std]-newStd) < 0:
			newStd = abs(0-errors[std])
		#correctedBelowStd.append(newStd)
		correctedBelowStd.append(errors[std] - belowStd[std])
	correctedAboveStd = []
	for std in range(0, len(aboveStd)):
		newStd = aboveStd[std]
		if errors[std]+newStd > 1 and labels[0] != 'Trees':
			newStd = abs(1-errors[std])
		#correctedAboveStd.append(newStd)
		correctedAboveStd.append(aboveStd[std] - errors[std])
	
		
	print correctedBelowStd
	print correctedAboveStd
	#Plot the error for the simulations
	xPositions = range(0, len(snpNums))
	p = ax.errorbar(xPositions, errors, yerr=[correctedBelowStd, correctedAboveStd], label='Normal', color=colors[colInd], linewidth=2)
	legendLines.append(p[0])
	
	#Repeat for the permutation run
	correctedBelowStdP = []
	for std in range(0, len(belowStdP)):

		correctedBelowStdP.append(pErrors[std] - belowStdP[std])
	correctedAboveStdP = []
	for std in range(0, len(aboveStdP)):
		
		#correctedAboveStd.append(newStd)
		correctedAboveStdP.append(aboveStdP[std] - pErrors[std])
		
	#Plot the error for the simulations
	xPositions = range(0, len(snpNums))
	print xPositions
	p = ax.errorbar(xPositions, pErrors, yerr=[correctedBelowStdP, correctedAboveStdP], label='Shuffled LAF', color=colors[colInd+1], linewidth=2)
	legendLines.append(p[0])
	
	
	#ax.set_ylim(lim[0],lim[1])
	ax.set_xlabel('Number of SNPs')
 	ax.set_ylabel('Error')
	#ax.set_xlim(-0.005,0.105)
	#extraticks = [0.005, 0.01, 0.015, 0.025, 0.03] #these are the ticks that are usually ignored
	#lim = ax.get_xlim()
	#print sorted(list(ax.get_xticks()) + extraticks)
	#ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	ax.set_xticks(xPositions)
	ax.set_xticklabels(snpNums, rotation=90)
	#ax.set_xlim(lim)
	plt.tight_layout()
	plt.legend()
	#plt.show()
	plt.savefig(title + '.svg')
	
	
	return 0



plotHorizontalDependencyInfluence(averagedCErrors, averagedPCErrors, groupedAboveStdC, groupedBelowStdC, groupedAboveStdCP, groupedBelowStdCP, snpNums, 'Copy numbers', 'Copy_numbers_hp', 0)
plotHorizontalDependencyInfluence(averagedAErrors, averagedPAErrors, groupedAboveStdA, groupedBelowStdA, groupedAboveStdAP, groupedBelowStdAP, snpNums, 'Alleles', 'Alleles_hp', 2)
plotHorizontalDependencyInfluence(averagedMuErrors, averagedPMuErrors, groupedAboveStdMu, groupedBelowStdMu, groupedAboveStdMuP, groupedBelowStdMuP, snpNums, 'Mu', 'Mu_hp', 4)
plotHorizontalDependencyInfluence(averagedTreeErrors, averagedPTreeErrors, groupedAboveStdT, groupedBelowStdT, groupedAboveStdTP, groupedBelowStdTP, snpNums, 'Trees', 'Trees_hp', 6)
	
plotHorizontalDependencyInfluence(averagedAncestrySwapErrors, averagedPAncestrySwapErrors, groupedAboveStdAncestry, groupedBelowStdAncestry, groupedAboveStdAncestryP, groupedBelowStdAncestry, snpNums, 'Trees', 'Ancestry_hp', 6)

#4. Another plot that we can make is show how much performance is gained with the horizontal dependency on the y axis, and the average distance between SNPs on the x axis (all in 1 plot is possible)



