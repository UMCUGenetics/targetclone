"""
	Script to visualize the effect of the influence of increasing numbers of SNVs on the performance of TargetClone. 

"""


import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

import computeTreeErrorOtherMetrics
from tree import Graph
import scipy as sp
import scipy.stats

#1. Define the folders containing the data that we want to visualize
motherFolder = sys.argv[1]
snvNum = [10, 50, 100, 500, 1000, 5000, 10000]

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
	
	groupedEuclideanErrors = dict()
	groupedAverageAncestrySwapErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + '/snvs_' + str(noiseLevel)
		
		#Read all the errors into one list for this noise level
		cErrors = []
		aErrors = []
		muErrors = []
		treeErrors = []
		ambiguityErrors = []
		ambiguityCorrectedErrors = []
		averagedAncestrySwapError = []
		
		treeSizes = []
		
		for subdir, dirs, files in os.walk(simulationFolder):
			
			
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
				if re.match('RealTrees', file): #read the file and obtain the error
					stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
					tree = eval(stringDict)
					realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
					treeSizes.append(len(realTree.edgeList))
				
				if re.match('EstimatedTrees', file): #read the file and obtain the error
					stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
					tree = eval(stringDict)
					inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
			
			[ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)
			
			#Instead of reporting the actual errors, what if we report percentages of how bad we could have done? 
			summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
			averagedAncestrySwapError.append(summedError / float(noOfSamplePairs))
					
			
		#Gather the data per noise level
		groupedCErrors[noiseLevel] = cErrors
		groupedAErrors[noiseLevel] = aErrors
		groupedMuErrors[noiseLevel] = muErrors
		groupedTreeErrors[noiseLevel] = treeErrors
		groupedAmbiguityErrors[noiseLevel] = ambiguityCorrectedErrors
		groupedAverageAncestrySwapErrors[noiseLevel] = averagedAncestrySwapError
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAverageAncestrySwapErrors]

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
		print noiseValues
		[mean, lower, upper] = mean_confidence_interval(noiseValues)
		
		groupedAboveStd.append(upper)
		groupedBelowStd.append(lower)
		
	sortedKeys, sortedBelow = zip(*sorted(zip(groupedErrors.keys(), groupedBelowStd)))
	sortedKeys, sortedAbove = zip(*sorted(zip(groupedErrors.keys(), groupedAboveStd)))
	return [sortedAbove, sortedBelow]

#Get the raw errors for the horizontal shuffle and the normal case
[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestrySwapErrors] = readDataIncludingPermutations(motherFolder, snvNum)

#Compute an average of the errors
print "C"
averagedCErrors = averageData(groupedCErrors, 'C')
print averagedCErrors
print "A"
averagedAErrors = averageData(groupedAErrors, 'C')
print averagedAErrors
print "Mu"
averagedMuErrors = averageData(groupedMuErrors, 'C')
print averagedMuErrors
print "T"
averagedTreeErrors = averageData(groupedTreeErrors, 'T')
print averagedTreeErrors
averagedAncestryErrors = averageData(groupedAncestrySwapErrors, 'C')
print averagedAncestryErrors


#Compute the standard deviation of the error (add later)
[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)
[groupedAboveStdAncestry, groupedBelowStdAncestry] = obtainStandardDeviations(groupedAncestrySwapErrors, averagedAncestryErrors)

#3. Make a plot with the error on the y axis, and the number of SNPs on the x axis. (4 plots per data type that we infer)

def plotHorizontalDependencyInfluence(errors, aboveStd, belowStd, snpNums, plotType, title, colInd, lim):
	
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
		
		
		correctedBelowStd.append(errors[std] - belowStd[std])
	correctedAboveStd = []
	for std in range(0, len(aboveStd)):
		
		correctedAboveStd.append(aboveStd[std] - errors[std])
		
	print correctedBelowStd
	print correctedAboveStd
	#Plot the error for the simulations
	xPositions = range(0, len(snpNums))
	p = ax.errorbar(xPositions, errors, yerr=[correctedBelowStd, correctedAboveStd], label='Normal', color=colors[colInd], linewidth=2)
	legendLines.append(p[0])
	
	ax.set_ylim(lim[0],lim[1])
	ax.set_xlabel('Number of SNVs')
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
	#plt.legend()
	#plt.show()
	plt.savefig(title + '.svg')
	
	
	return 0




plotHorizontalDependencyInfluence(averagedCErrors, groupedAboveStdC, groupedBelowStdC, snvNum, 'Copy numbers', 'Copy_numbers_snv', 0, [-0.05,0.4])
plotHorizontalDependencyInfluence(averagedAErrors, groupedAboveStdA, groupedBelowStdA, snvNum, 'Alleles', 'Alleles_snv', 2, [-0.05,0.4])
plotHorizontalDependencyInfluence(averagedMuErrors, groupedAboveStdMu, groupedBelowStdMu, snvNum, 'Mu', 'Mu_snv', 4, [-0.05,0.3])
plotHorizontalDependencyInfluence(averagedTreeErrors, groupedAboveStdT, groupedBelowStdT, snvNum, 'Trees', 'Trees_snv', 6, [-0.05,6])
plotHorizontalDependencyInfluence(averagedAncestryErrors, groupedAboveStdAncestry, groupedBelowStdAncestry, snvNum, 'Trees', 'Ancestry_snv', 6, [0,0.3])


