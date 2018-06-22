"""
	Script to visualize the effect of the influence of increasing numbers of SNVs on the performance of TargetClone. 

"""


import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

#1. Define the folders containing the data that we want to visualize
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
	
	groupedEuclideanErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + '/snvs' + str(noiseLevel)
		
		#Read all the errors into one list for this noise level
		cErrors = []
		aErrors = []
		muErrors = []
		treeErrors = []
		ambiguityErrors = []
		ambiguityCorrectedErrors = []
		
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
						
			
		#Gather the data per noise level
		groupedCErrors[noiseLevel] = cErrors
		groupedAErrors[noiseLevel] = aErrors
		groupedMuErrors[noiseLevel] = muErrors
		groupedTreeErrors[noiseLevel] = treeErrors
		groupedAmbiguityErrors[noiseLevel] = ambiguityCorrectedErrors
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors]

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

def obtainStandardDeviations(groupedErrors, averagedError):
	
	#Compute standard deviations, see if it is better
	groupedAboveStd = []
	groupedBelowStd = []
	
	#It is more interesting to show the quantiles rather than the standard deviation
	for noiseLevelInd in range(0, len(groupedErrors.keys())):
		noiseValues = groupedErrors[groupedErrors.keys()[noiseLevelInd]]
		q1 = np.std(noiseValues)
		q3 = np.std(noiseValues)
		groupedAboveStd.append(q3)
		groupedBelowStd.append(q1)
		
	sortedKeys, sortedBelow = zip(*sorted(zip(groupedErrors.keys(), groupedBelowStd)))
	sortedKeys, sortedAbove = zip(*sorted(zip(groupedErrors.keys(), groupedAboveStd)))
	return [sortedAbove, sortedBelow]

#Get the raw errors for the horizontal shuffle and the normal case
[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors] = readDataIncludingPermutations(motherFolder, snpNums)

#Compute an average of the errors
print "C"
averagedCErrors = averageData(groupedCErrors, 'C')
print "A"
averagedAErrors = averageData(groupedAErrors, 'C')
print "Mu"
averagedMuErrors = averageData(groupedMuErrors, 'C')
print "T"
averagedTreeErrors = averageData(groupedTreeErrors, 'T')
print averagedTreeErrors


#Compute the standard deviation of the error (add later)
[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)


#3. Make a plot with the error on the y axis, and the number of SNPs on the x axis. (4 plots per data type that we infer)

def plotHorizontalDependencyInfluence(errors, aboveStd, belowStd, snpNums, plotType, title, colInd):
	
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
		correctedBelowStd.append(newStd)
	correctedAboveStd = []
	for std in range(0, len(aboveStd)):
		newStd = aboveStd[std]
		if errors[std]+newStd > 1 and plotType != 'Trees':
			newStd = abs(1-errors[std])
		correctedAboveStd.append(newStd)
	
		
	print correctedBelowStd
	print correctedAboveStd
	#Plot the error for the simulations
	xPositions = range(0, len(snpNums))
	p = ax.errorbar(xPositions, errors, yerr=[correctedBelowStd, correctedAboveStd], label='Normal', color=colors[colInd], linewidth=2)
	legendLines.append(p[0])
	
	#ax.set_ylim(lim[0],lim[1])
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
	plt.legend()
	plt.show()
	#plt.savefig(title + '.svg')
	
	
	return 0


snvNum = [10, 50, 100, 500, 1000, 5000, 10000, 50000]

plotHorizontalDependencyInfluence(averagedCErrors, groupedAboveStdC, groupedBelowStdC, snpNums, 'Copy numbers', 'Copy_numbers_hp', 0)
plotHorizontalDependencyInfluence(averagedAErrors, groupedAboveStdA, groupedBelowStdA, snpNums, 'Alleles', 'Alleles_hp', 2)
plotHorizontalDependencyInfluence(averagedMuErrors, groupedAboveStdMu, groupedBelowStdMu, snpNums, 'Mu', 'Mu_hp', 4)
plotHorizontalDependencyInfluence(averagedTreeErrors, groupedAboveStdT, groupedBelowStdT, snpNums, 'Trees', 'Trees_hp', 6)
	



