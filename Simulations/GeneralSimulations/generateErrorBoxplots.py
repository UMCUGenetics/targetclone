import sys
sys.path.insert(0, '../../TargetClone/')

import os
import re
import ast
import matplotlib.pyplot as plt
from pylab import axes
import numpy as np
from distance import EventDistances
import computeTreeErrorOtherMetrics
from tree import Graph
from combinations import AlleleCombination
from alleles import Alleles
from mu import Mu
from math import sqrt
from copy import deepcopy
from scipy import stats

import pickle as pkl

plt.rcParams.update({'font.size': 14})




def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines

def readDataIncludingPermutations(dataFolder, noiseLevels, addition):
	
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
	groupedAncestrySwapErrorsAbsentInInferred = dict()
	groupedAncestrySwapErrorsPresentInInferred = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + noiseLevel
		
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
		
		euclideanErrors = []
		ancestrySwapErrorsAbsentInInferred = []
		ancestrySwapErrorsPresentInInferred = []
		
		realTree = None
		inferredTree = None
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

				#Repeat for the permutation errors
				if re.match('pCError', file): #read the file and obtain the error
					pCErrors += collectErrorsFromFile(file, subdir)
				if re.match('pAError', file): #read the file and obtain the error
					pAErrors += collectErrorsFromFile(file, subdir)
				if re.match('pMuError', file): #read the file and obtain the error
					pMuErrors += collectErrorsFromFile(file, subdir)
				if re.match('pTreeError', file): #read the file and obtain the error
					pTreeErrors += collectErrorsFromFile(file, subdir)
				
				#Also collect the real tree and inferred tree to compute the anscestry swap errors
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
			
			#Instead of reporting the actual errors, what if we report percentages of how bad we could have done? 
			
			ancestrySwapErrorsAbsentInInferred.append(ancestrySwapErrorAbsentInInferred / float(noOfSamplePairs))
			ancestrySwapErrorsPresentInInferred.append(ancestrySwapErrorPresentInInferred / float(noOfSamplePairs))
		
		print "tree sizes:"
		print sum(treeSizes) / float(len(treeSizes))
		#exit()
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
		groupedAncestrySwapErrorsAbsentInInferred[noiseLevel] = ancestrySwapErrorsAbsentInInferred
		groupedAncestrySwapErrorsPresentInInferred[noiseLevel] = ancestrySwapErrorsPresentInInferred
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, groupedAncestrySwapErrorsAbsentInInferred, groupedAncestrySwapErrorsPresentInInferred]

def obtainStandardDeviations(groupedErrors, averagedError):
	
	
	#Make sure that we sort the keys correctly
	
	# groupedAboveStd = []
	# groupedBelowStd = []
	# 
	# #It is more interesting to show the quantiles rather than the standard deviation
	# for noiseLevelInd in range(0, len(groupedErrors.keys())):
	# 	noiseValues = groupedErrors[groupedErrors.keys()[noiseLevelInd]]
	# 	q1 = np.percentile(noiseValues, 25)
	# 	q3 = np.percentile(noiseValues, 75)
	# 	groupedAboveStd.append(q3)
	# 	groupedBelowStd.append(q1)
	# 	
	# sortedKeys, sortedBelow = zip(*sorted(zip(groupedErrors.keys(), groupedBelowStd)))
	# sortedKeys, sortedAbove = zip(*sorted(zip(groupedErrors.keys(), groupedAboveStd)))
	# return [sortedAbove, sortedBelow]


	#Compute standard deviations, see if it is better
	groupedAboveStd = []
	groupedBelowStd = []
	
	#It is more interesting to show the quantiles rather than the standard deviation
	for noiseLevelInd in range(0, len(groupedErrors.keys())):
		noiseValues = groupedErrors[groupedErrors.keys()[noiseLevelInd]]
		currentStd = np.std(noiseValues)
		currentMean = np.mean(noiseValues)
		conf_int = stats.norm.interval(0.95, loc=currentMean, scale=currentStd)
		print "noise level: ", groupedErrors.keys()[noiseLevelInd]
		print "current mean: ", currentMean
		print "current std: ", currentStd
		print "confidence interval: ", conf_int
		#q1 = np.percentile(noiseValues, 5)
		#q3 = np.percentile(noiseValues, 95)
		#p025 = df.groupby('category')['number'].quantile(0.025)
		#p975 = df.groupby('category')['number'].quantile(0.975)
		
		groupedAboveStd.append(conf_int[1])
		groupedBelowStd.append(conf_int[0])
		
	sortedKeys, sortedBelow = zip(*sorted(zip(groupedErrors.keys(), groupedBelowStd)))
	sortedKeys, sortedAbove = zip(*sorted(zip(groupedErrors.keys(), groupedAboveStd)))
	
	return [sortedAbove, sortedBelow]

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


def averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, swapErrorAbsent, swapErrorPresent):

	averagedCErrors = averageData(groupedCErrors, 'C')
	averagedAErrors = averageData(groupedAErrors, 'A')
	averagedMuErrors = averageData(groupedMuErrors, 'Mu')
	
	averagedTreeErrors = averageData(groupedTreeErrors, 'T')
	averagedAmbiguityErrors = averageData(groupedAmbiguityErrors, 'C')

	averagedSwapErrorsAbsent = averageData(swapErrorAbsent, 'T')
	averagedSwapErrorsPresent = averageData(swapErrorPresent, 'T')

	
	return [averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors, averagedSwapErrorsAbsent, averagedSwapErrorsPresent]

def computeRandomCaseError(dataFolder):
	
	#F1. Read the data from the data folder
	#Make dummy noise levels to read the right folder
	noiseLevels = ['generic_random']
	[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, ancestrySwapErrorAbsent, ancestrySwapErrorPresent] = \
	readDataIncludingPermutations(dataFolder, noiseLevels, '')
	
	
	print "distribution of random absent errors: ", ancestrySwapErrorAbsent
	print "distribution of random present errors: ", ancestrySwapErrorPresent
	
	#make a plot to see how the errors are distributed
	#The distribution of the z-score is the same as the distribution of the permutations! It is only standardized. 

	#F2. Average the errors (one noise level)
	[averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, a, averagedSwapErrorsAbsent, averagedSwapErrorsPresent] = \
	averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedCErrors, ancestrySwapErrorAbsent, ancestrySwapErrorPresent)
	
	#Obtain the standard deviations
	
	print "C std: "
	randomCStd = obtainStandardDeviations(groupedCErrors, averagedCErrors)
	print "A std: "
	randomAStd = obtainStandardDeviations(groupedAErrors, averagedAErrors)
	print "Mu std: "
	randomMuStd = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
	randomTreeStd = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)
	swapAbsentStd = obtainStandardDeviations(ancestrySwapErrorAbsent, averagedSwapErrorsAbsent)
	swapPresentStd = obtainStandardDeviations(ancestrySwapErrorPresent, averagedSwapErrorsPresent)
	
	#randomCStd = np.std(groupedCErrors[noiseLevels[0]])
	#randomAStd = np.std(groupedAErrors[noiseLevels[0]])
	#randomMuStd = np.std(groupedMuErrors[noiseLevels[0]])
	#randomTreeStd = np.std(groupedTreeErrors[noiseLevels[0]])

	#F3. Print the errors.
	print averagedCErrors
	print averagedAErrors
	print averagedMuErrors
	print averagedTreeErrors
	
	print averagedSwapErrorsAbsent
	print averagedSwapErrorsPresent
	
	return [averagedCErrors[0], averagedAErrors[0], averagedMuErrors[0], averagedTreeErrors[0], averagedSwapErrorsAbsent[0], averagedSwapErrorsPresent[0], randomCStd, randomAStd, randomMuStd, randomTreeStd, swapAbsentStd, swapPresentStd]

#2. Calculate statistics for the completely random case
dataFolder = 'Results/'
[randomCError, randomAError, randomMuError, randomTreeError, randomSwapErrorAbsent, randomSwapErrorPresent, randomCStd, randomAStd, randomMuStd, randomTreeStd, swapAbsentStd, swapPresentStd] = computeRandomCaseError(dataFolder)


def readData(dataFolder, noiseLevels, addition):
	
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
	groupedAncestryAbsentErrors = dict()
	groupedAncestryPresentErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + str(noiseLevel) + addition
		
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
		
		euclideanErrors = []
		ancestrySwapErrorsAbsentInInferred = []
		ancestrySwapErrorsPresentInInferred = []
		
		lafMatrix = None
		realTree = None
		snvMatrix = None
		
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
				if re.match('ambiguityError', file):
					ambiguityErrors += collectErrorsFromFile(file, subdir)
				if re.match('ambiguityCorrectedError', file):
					ambiguityCorrectedErrors += collectErrorsFromFile(file, subdir)
					
					
				#Repeat for the permutation errors
				if re.match('pCError', file): #read the file and obtain the error
					pCErrors += collectErrorsFromFile(file, subdir)
				if re.match('pAError', file): #read the file and obtain the error
					pAErrors += collectErrorsFromFile(file, subdir)
				if re.match('pMuError', file): #read the file and obtain the error
					pMuErrors += collectErrorsFromFile(file, subdir)
				if re.match('pTreeError', file): #read the file and obtain the error
					pTreeErrors += collectErrorsFromFile(file, subdir)
				if re.match('pAmbiguityError', file):
					pAmbiguityErrors += collectErrorsFromFile(file, subdir)
				if re.match('pAmbiguityCorrectedError', file):
					pAmbiguityCorrectedErrors += collectErrorsFromFile(file, subdir)
					
				#Also collect the real tree and inferred tree to compute the anscestry swap errors
				if re.match('RealTrees', file): #read the file and obtain the error
					stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
					tree = eval(stringDict)
					realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
				
				if re.match('EstimatedTrees', file): #read the file and obtain the error
					stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
					tree = eval(stringDict)
					inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
			
			
			#Compute the ancestry swap error
			[ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)
			
			ancestrySwapErrorsAbsentInInferred.append(ancestrySwapErrorAbsentInInferred / float(noOfSamplePairs))
			ancestrySwapErrorsPresentInInferred.append(ancestrySwapErrorPresentInInferred / float(noOfSamplePairs))
				
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
		groupedAncestryAbsentErrors[noiseLevel] = ancestrySwapErrorsAbsentInInferred
		groupedAncestryPresentErrors[noiseLevel] = ancestrySwapErrorsPresentInInferred
		
		#Move this to a function to make it better
		#Also compute the Euclidean distance trees for each noise levels, add this as an additional error
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors]

def convertToZScores(randomError, randomStd, realErrors):
	
	standardizedErrors = dict()
	
	for noiseLevel in realErrors.keys():
		
		values = realErrors[noiseLevel]
		print randomError
		print randomStd
		print values
		z = [(realError - randomError) / randomStd for realError in values]
		standardizedErrors[noiseLevel] = z 
	
	return standardizedErrors

def plotData(noiseLevels, errors, aboveStd, belowStd, randomError, randomStd, labels, colInd, lim, title):
	
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	plt.figure()
	ax = plt.gca()
	
	# 
	# 
	# stdAbove = [randomStd[0][0] for i in [randomError]*len(noiseLevels)]
	# stdBelow = [randomStd[1][0] for i in [randomError]*len(noiseLevels)]
	# print "random error stds: "
	# print stdAbove
	# print stdBelow
	# 
	# 
	# ax.plot(noiseLevels, [randomError]*len(noiseLevels), color="#CDCDCD")
	# ax.plot(noiseLevels, stdAbove, color="#CDCDCD")
	# ax.plot(noiseLevels, stdBelow, color="#CDCDCD")
	# 
	# ax.fill_between(noiseLevels, stdBelow, stdAbove, alpha=0.5, edgecolor='#CDCDCD', facecolor='#CDCDCD')
	# 
	# legendLines = []
	
	correctedBelowStd = []
	for std in range(0, len(belowStd)):
		newStd = belowStd[std]
		if (errors[std]-newStd) < 0:
			newStd = abs(0-errors[std])
		#correctedBelowStd.append(newStd)
		correctedBelowStd.append(belowStd[std])
	correctedAboveStd = []
	for std in range(0, len(aboveStd)):
		newStd = aboveStd[std]
		if errors[std]+newStd > 1 and labels[0] != 'Trees':
			newStd = abs(1-errors[std])
		#correctedAboveStd.append(newStd)
		correctedAboveStd.append(aboveStd[std])
	
	print "corrected stds: "	
	print correctedBelowStd
	print correctedAboveStd
	#Plot the error for the simulations
	
	print "plotting means: ", errors
	print "plotting above: ", aboveStd
	print "plotting below: ", belowStd
	
	p = ax.errorbar(noiseLevels, errors, yerr=[belowStd, aboveStd], label='$E_C$', color=colors[colInd], linewidth=2)
	#legendLines.append(p[0])
	
	#ax.legend(legendLines, labels, loc=2, numpoints=1)
	
	ax.set_ylim(lim[0],lim[1])
	ax.set_xlabel(r'Sequencing noise ($\sigma$)')
 	ax.set_ylabel('Error')
	ax.set_xlim(-0.005,0.105)
	extraticks = [0.005, 0.01, 0.015, 0.025, 0.03] #these are the ticks that are usually ignored
	lim = ax.get_xlim()
	#print sorted(list(ax.get_xticks()) + extraticks)
	ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	ax.set_xticklabels(['', 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1, ''], rotation=90)
	ax.set_xlim(lim)
	plt.tight_layout()
	
	#plt.show()
	plt.savefig(title)


def plotFigureOne(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors, averagedSwapAbsentErrors, averagedSwapPresentErrors, 
				  noiseLevels, groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu,
				  groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT, ambiguityScores, groupedAboveStdAmb, groupedBelowStdAmb,
				  ambiguityScoresRandom, groupedAboveStdAmbR, groupedBelowStdAmbR,
				  groupedAboveStdSwapAbsent, groupedBelowStdSwapAbsent, groupedAboveStdSwapPresent, groupedBelowStdSwapPresent):
	
	
	
	ambiguityErrors = []
	ambiguityErrorsRandom = 1 - ambiguityScoresRandom
	for score in ambiguityScores:
		ambiguityErrors.append(1-score)
	
	#Also swap the ambiguity stds?
	plotData(noiseLevels, averagedCErrors, groupedAboveStdC, groupedBelowStdC, randomCError, randomCStd, ['Copy numbers'], 0, [0,1], 'fig3_C.svg')
	plotData(noiseLevels, averagedAErrors, groupedAboveStdA, groupedBelowStdA, randomAError, randomAStd, ['Alleles'], 1, [0,1], 'fig3_A.svg')
	plotData(noiseLevels, averagedMuErrors, groupedAboveStdMu, groupedBelowStdMu, randomMuError, randomMuStd, ['Tumor fraction'], 3, [0,0.6], 'fig3_Mu.svg')
	plotData(noiseLevels, averagedTreeErrors, groupedAboveStdT, groupedBelowStdT, randomTreeError, randomTreeStd, ['Trees'], 4, [-1,10], 'fig3_T.svg')
	plotData(noiseLevels, averagedAmbiguityErrors, groupedAboveStdAmb, groupedBelowStdAmb, ambiguityErrorsRandom, [groupedAboveStdAmbR,groupedBelowStdAmbR], ['Resolved ambiguities'], 5, [0,1], 'fig3_Amb.svg')
	
	#Swap errors.Should these be visualized in the same figure, or keep them separate?
	plotData(noiseLevels, averagedSwapAbsentErrors, groupedAboveStdSwapAbsent, groupedBelowStdSwapAbsent, randomSwapErrorAbsent, swapAbsentStd, ['Ancestry swap'], 5, [0,1], 'fig3_AncestrySwapAbsent.svg')
	plotData(noiseLevels, averagedSwapPresentErrors, groupedAboveStdSwapPresent, groupedBelowStdSwapPresent, randomSwapErrorPresent, swapPresentStd, ['Ancestry swap'], 6, [0,1], 'fig3_AncestrySwapPresent.svg')

def plotBoxplots(noiseLevels, errors, randomErrors, title):
	
	randomErrorList = randomErrors['generic_random']
	print randomErrorList
	
	#1. Make a list where each entry is the errors at a different noise level
	#Sort the error list by their noise value order
	
	errorList = []
	for noiseLevel in noiseLevels:
		errorList.append(errors[noiseLevel])
	
	errorList.append(randomErrorList)
	
	
	ax = plt.gca()
	ax.boxplot(errorList)
	
	
	xtickLabels = deepcopy(noiseLevels)
	xtickLabels.append('random')
	ax.set_xticklabels(xtickLabels, rotation=90)
	
	plt.savefig(title)
	plt.show()
	

def makeBoxPlotFigure(dataFolder, noiseLevels):
	
	#0. Obtain the raw error data for the random case
	randomNoiseLevels = ['generic_random']
	randomDataFolder = 'Results/'
	[randomCErrors, randomAErrors, randomMuErrors, randomTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, ancestrySwapErrorAbsentRandom, ancestrySwapErrorPresentRandom] = \
	readDataIncludingPermutations(randomDataFolder, randomNoiseLevels, '')
	
	#F1. Read the data from all the simulation folders (The normal and permuted errors)
	[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors] = readData(dataFolder, noiseLevels, '')
	
	plotBoxplots(noiseLevels, groupedCErrors, randomCErrors, 'boxplot_C.svg')
	plotBoxplots(noiseLevels, groupedAErrors, randomAErrors, 'boxplot_A.svg')
	plotBoxplots(noiseLevels, groupedMuErrors, randomMuErrors, 'boxplot_Mu.svg')
	plotBoxplots(noiseLevels, groupedTreeErrors, randomTreeErrors, 'boxplot_T.svg')
	
	plotBoxplots(noiseLevels, groupedAncestryAbsentErrors, ancestrySwapErrorAbsentRandom, 'boxplot_AncestrySwapAbsent.svg')
	plotBoxplots(noiseLevels, groupedAncestryPresentErrors, ancestrySwapErrorPresentRandom, 'boxplot_AncestrySwapPresent.svg')
	
	return 0

	


def plotFigureOneOld(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors,
				  noiseLevels, groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu,
				  groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT, ambiguityScores, groupedAboveStdAmb, groupedBelowStdAmb,
				  ambiguityScoresRandom, groupedAboveStdAmbR, groupedBelowStdAmbR):
	
	#Make a plot of the noise levels vs the error
	#Add all lines in one figure, but keep the permutation one separate (?)
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	
	ax = plt.gca()
	ax2 = ax.twinx()
	
	#Convert the ambiguityScores to errors
	ambiguityErrors = []
	ambiguityErrorsRandom = []
	for score in ambiguityScores:
		ambiguityErrors.append(1-score)
	for score in ambiguityScoresRandom:
		ambiguityErrorsRandom.append(1-score)
		
		
	#Make figures one by one in separate plots
	
	#1. Add a line for the C error
	
	legendLines = []
	legendLines2 = []
	labels = ['Copy numbers', 'Alleles', 'Tumor fraction', 'Resolved ambiguities']
	#Plot the error for the simulations
	p = ax.errorbar(noiseLevels, averagedCErrors, yerr=[groupedBelowStdC, groupedAboveStdC], label='$E_C$', color=colors[0], linewidth=2)
	legendLines.append(p[0])
	
	
	
	#Add some offsets to the values
	aOffset = 0.0003
	muOffset = 0.0007
	treeOffset = 0.0011
	ambOffset = 0.00015
	
	noiseLevelsA = [i-aOffset for i in noiseLevels]
	noiseLevelsMu = [i-muOffset for i in noiseLevels]
	noiseLevelsTree = [i-treeOffset for i in noiseLevels]
	noiseLevelsAmb = [i-ambOffset for i in noiseLevels]
	
	#Add the lines for the permutations (first try with one)
	stdAbove = [i + randomCStd[0][0] for i in [randomCError]*len(noiseLevels)]

	stdBelow = [i - randomCStd[1][0] for i in [randomCError]*len(noiseLevels)]
	
	ax.plot(noiseLevels, [randomCError]*len(noiseLevels), color="#CC4F1B")
	ax.plot(noiseLevels, stdAbove, color="#CC4F1B")
	ax.plot(noiseLevels, stdBelow, color="#CC4F1B")
	
	ax.fill_between(noiseLevels, stdBelow, stdAbove, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
	
	#1. Add a line for the C error
	
	legendLines = []
	legendLines2 = []
	labels = ['Copy numbers', 'Alleles', 'Tumor fraction', 'Resolved ambiguities']
	#Plot the error for the simulations
	p = ax.errorbar(noiseLevels, averagedCErrors, yerr=[groupedBelowStdC, groupedAboveStdC], label='$E_C$', color=colors[0], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevelsA, averagedAErrors, yerr=[groupedBelowStdA, groupedAboveStdA], label='$E_A$', color=colors[1], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevelsMu, averagedMuErrors, yerr=[groupedBelowStdMu, groupedAboveStdMu], label='r$E_\mu$', color=colors[3], linewidth=2)
	legendLines.append(p[0])
	p = ax2.errorbar(noiseLevelsTree, averagedTreeErrors, yerr=[groupedBelowStdT, groupedAboveStdT], label='$E_T$', color=colors[4], linewidth=2)
	legendLines2.append(p[0])

	p = ax.errorbar(noiseLevelsAmb, ambiguityErrors, yerr=[groupedBelowStdAmb, groupedAboveStdAmb], label='r$E_\mu$', color=colors[5], linewidth=2)
	legendLines.append(p[0])
	
	labels2 = ['Trees']
	ax.legend(legendLines, labels, loc=2, numpoints=1)
	ax2.legend(legendLines2, labels2)
	ax2.set_ylim(0,6)
	ax.set_ylim(0,1.5)
	ax.set_xlabel(r'Sequencing noise ($\sigma$)')
 	ax.set_ylabel('Error')
	ax2.set_ylabel('Error')
	ax.set_xlim(-0.005,0.105)
	extraticks = [0.005, 0.01, 0.015, 0.025, 0.03] #these are the ticks that are usually ignored
	lim = ax.get_xlim()
	#print sorted(list(ax.get_xticks()) + extraticks)
	ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	ax.set_xticklabels(['', 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1, ''], rotation=90)
	ax.set_xlim(lim)
	plt.tight_layout()
	
	plt.show()
	#plt.savefig('errorFunctionOfNoise.svg')

#Compute how many positions are actually ambiguous

def computeAmbiguousPositions():
	
	#We can generate the allele list with the event distances function
	kmin = 1 #the kmin and kmax used in the simulations
	kmax = 6 
	eventDistance = EventDistances(kmin, kmax)
	
	#get the allele list
	alleleList = eventDistance.alleleList
	
	#make sure that alleles are not duplicated
	LAFAndCombinations = dict()
	normalAlleles = Alleles(1,1)
	for allele in alleleList:
		AOccurrences = [m.start() for m in re.finditer('A', allele)]
		ACount = len(AOccurrences)
		BOccurrences = [m.start() for m in re.finditer('B', allele)]
		BCount = len(BOccurrences)
		
		alleleObj = Alleles(ACount, BCount)
		if BCount > ACount or BCount == ACount:
			alleleCombination = AlleleCombination([normalAlleles, alleleObj])
			
			for muIndex in range(0,101):
				
				LAF = alleleCombination.computeLAF(Mu(muIndex)) #compute the LAF that each combination would generate
				if LAF not in LAFAndCombinations.keys():
					#LAFAndCombinations[LAF] = []
					LAFAndCombinations[LAF] = 0
				
				#LAFAndCombinations[LAF].append((alleleObj.getAllelesAsString(), muIndex))
				LAFAndCombinations[LAF] += 1
			
	
	#print LAFAndCombinations
	#For every mu, we should check which LAF the combination with normal would generate
	
	#With this dictionary, we can check if a LAF has more than one solution. If true, then we can check and see if the position is correct or not. With that, we compute a score showing the
	#number of ambiguous positions that we were able to infer correctly. This score is a higher wow factor than 
	return LAFAndCombinations
	
	
	
LAFAndCombinations = computeAmbiguousPositions() #This is a pre-computed dictionary of all LAF and the associated number of combinations.

#For each simulated dataset, we look at the real mu and the real A and compute the LAF. Does it correspond to a LAF that has multiple solutions?
#If yes, did we get the position correct?
#Perhaps, if no, what is the reason? Is it because we inferred another ambiguous solution, or did the method go wrong completely? (This can maybe be interesting as well)

def computeCorrectAmbiguityScore(LAFAndCombinations, simulationFolder):
	ambiguityScores = []
	ambiguities = []
	correctAmbiguityPositions = 0
	totalAmbiguousPositions = 0
	totalSize = 0
	#We need to read the actual A matrix values and also the mu
	normalAlleles = Alleles(1,1)
	#1. read the simulated A matrix
	for subdir, dirs, files in os.walk(simulationFolder):
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		for file in files:
			if re.match('RealA', file): #read the file and obtain the a matrix
				realAMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
			if re.match('RealMu', file): #also read the real mu
				realMu = collectErrorsFromFile(file, subdir)
				
			#Then load the inferred A and mu
			if re.match('EstimatedA', file): #read the file and obtain the a matrix
				estimatedAMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
			if re.match('EstimatedMu', file): #also read the real mu
				estimatedMu = collectErrorsFromFile(file, subdir)
			
		#Compute the LAF that each measurement in the real data would generate
		for row in range(0, realAMatrix.shape[0]):
			for col in range(0, realAMatrix.shape[1]):
				totalSize += 1
				#generate allele object
				allele = realAMatrix[row][col]
				AOccurrences = [m.start() for m in re.finditer('A', allele)]
				ACount = len(AOccurrences)
				BOccurrences = [m.start() for m in re.finditer('B', allele)]
				BCount = len(BOccurrences)
				
				alleleObj = Alleles(ACount, BCount)
				alleleCombination = AlleleCombination([normalAlleles, alleleObj])
				
				#Compute the LAF this combination would generate
				muNormal = 1-(realMu[col])
				realMuObj = Mu(int(muNormal*100)) #this function only takes integer indices!
				realLAF = alleleCombination.computeLAF(realMuObj)
				
				#Check if this LAF is ambiguous y/n.
				ambiguousCount = LAFAndCombinations[realLAF]
				
				#If the ambiguous count > 1 and we are correct, we make a note of that.
				if ambiguousCount > 1:
					totalAmbiguousPositions += 1
					if realAMatrix[row][col] == estimatedAMatrix[row][col]:
						correctAmbiguityPositions += 1
				
		#Divide the ambiguity score by the total number of positions.
		#print correctAmbiguityPositions
		#print correctAmbiguityPositions / float(totalAmbiguousPositions)
		#print totalAmbiguousPositions / float(totalSize)
		#ambiguityScores.append(correctAmbiguityPositions / float(totalAmbiguousPositions)) #Reporting as % of ambiuguities
		#Reporting the ambiguity scores as the fraction of the total
		ambiguityScores.append(correctAmbiguityPositions / float(totalSize))
		
		ambiguities.append(totalAmbiguousPositions / float(totalSize))
	
	#Compute an average for every noise level.
	
	#convert to z-scores
	
	averageAmbiguityScore = sum(ambiguityScores) / float(len(ambiguityScores))
	averageAmbiguities = sum(ambiguities) / float(len(ambiguities))
	
	return [averageAmbiguities, averageAmbiguityScore, ambiguityScores]

def generateFigureOne(dataFolder, noiseLevels, ambiguityScores, groupedAmbiguities, averageAmbiguityScoreRandom, groupedAmbiguityScoresRandom):
	
	#F1. Read the data from all the simulation folders (The normal and permuted errors)
	[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors] = readData(dataFolder, noiseLevels, '')
	# 
	print groupedCErrors
	# print np.percentile(groupedCErrors[0.01], 0)
	# print np.percentile(groupedCErrors[0.01], 25)
	# print np.percentile(groupedCErrors[0.01], 50)
	# print np.percentile(groupedCErrors[0.01], 75)
	# print np.percentile(groupedCErrors[0.01], 100)
	
	
	
	
	
	#F2. Average the errors per noise level
	[averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors, averagedSwapAbsentErrors, averagedSwapPresentErrors] = \
	averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors)
	
	
	#Obtain the standard deviation above and below the mean for the averages
	print "C std:"
	[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
	print "A std: "
	[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
	print "mu std: "
	[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
	[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)
	
	[groupedAboveStdAmb, groupedBelowStdAmb] = obtainStandardDeviations(groupedAmbiguities, ambiguityScores)

	[groupedAboveStdAmbR, groupedBelowStdAmbR] = obtainStandardDeviations(groupedAmbiguityScoresRandom, averageAmbiguityScoreRandom)
	
	[groupedAboveStdSwapAbsent, groupedBelowStdSwapAbsent] = obtainStandardDeviations(groupedAncestryAbsentErrors, averagedSwapAbsentErrors)
	[groupedAboveStdSwapPresent, groupedBelowStdSwapPresent] = obtainStandardDeviations(groupedAncestryPresentErrors, averagedSwapPresentErrors)

	#F3. Plot the error per noise level in one figure
	
	plotFigureOne(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors, averagedSwapAbsentErrors, averagedSwapPresentErrors, 
				  noiseLevels, groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu, groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT,
				  ambiguityScores, groupedAboveStdAmb, groupedBelowStdAmb, averageAmbiguityScoreRandom[0], groupedAboveStdAmbR, groupedBelowStdAmbR,
				  groupedAboveStdSwapAbsent, groupedBelowStdSwapAbsent, groupedAboveStdSwapPresent, groupedBelowStdSwapPresent)
	
	return 0


# #1. Make figure one panels A-D
# #For figure one, we average the errors across the noise levels and plot all values together in one tree.
# #For every simulation run, we also generate euclidean trees and compute the error. 
simulationFolder = 'Results/generic_noise'
dataFolder = 'Results/generic_noise'
#noiseLevels = [0, 0.005, 0.01]
noiseLevels = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1]
noiseLevels = [0.015, 0.02, 0.025]
#noiseLevels = [0]
#Put this here so that we can easily re-use it for other figures
ambiguities = []
ambiguityScores = []
groupedAmbiguities = dict()

for noiseLevel in noiseLevels:
	currentSimulationFolder = simulationFolder + str(noiseLevel) + '/'
	print "noise level: ", noiseLevel
	[averageAmbiguities, averageAmbiguityScore, allAmbiguityScores] = computeCorrectAmbiguityScore(LAFAndCombinations, currentSimulationFolder)
	groupedAmbiguities[noiseLevel] = allAmbiguityScores
	ambiguities.append(averageAmbiguities)
	ambiguityScores.append(averageAmbiguityScore)
	
#Also compute the ambiguity scores in the permuted data

permutationFolder = 'Results/generic_random/'

[averageAmbiguitiesRandom, averageAmbiguityScoreRandom, allAmbiguityScoresRandom] = computeCorrectAmbiguityScore(LAFAndCombinations, permutationFolder)
groupedAmbiguitiesRandom = dict()
groupedAmbiguitiesRandom[0] = allAmbiguityScoresRandom
generateFigureOne(dataFolder, noiseLevels, ambiguityScores, groupedAmbiguities, [averageAmbiguityScoreRandom], groupedAmbiguitiesRandom)
#makeBoxPlotFigure(dataFolder, noiseLevels)

exit()

#Make the random restarts figure


def plotDataRestarts(noiseLevels, errors, aboveStd, belowStd, labels, colInd, lim, title):
	
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
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
		if errors[std]+newStd > 1 and labels[0] != 'Trees':
			newStd = abs(1-errors[std])
		correctedAboveStd.append(newStd)
	
		
	print correctedBelowStd
	print correctedAboveStd
	#Plot the error for the simulations
	p = ax.errorbar(noiseLevels, errors, yerr=[correctedBelowStd, correctedAboveStd], label='$E_C$', color=colors[colInd], linewidth=2)
	legendLines.append(p[0])
	
	#ax.legend(legendLines, labels, loc=2, numpoints=1)
	
	ax.set_ylim(lim[0],lim[1])
	ax.set_xlabel(r'Sequencing noise ($\sigma$)')
 	ax.set_ylabel('Error')
	ax.set_xlim(0.00,0.035)
	#extraticks = [0.005, 0.01, 0.015, 0.025, 0.03] #these are the ticks that are usually ignored
	lim = ax.get_xlim()
	#print sorted(list(ax.get_xticks()) + extraticks)
	#ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	ax.set_xticklabels(['', 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, ''], rotation=90)
	ax.set_xlim(lim)
	plt.tight_layout()
	
	#plt.show()
	plt.savefig(title)
	
def plotFigureOneRandomRestarts(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors,
				  noiseLevels, groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu,
				  groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT):
	
	plotDataRestarts(noiseLevels, averagedCErrors, groupedAboveStdC, groupedBelowStdC, ['Copy numbers'], 0, [0,1], 'randomRestarts_C.svg')
	plotDataRestarts(noiseLevels, averagedAErrors, groupedAboveStdA, groupedBelowStdA, ['Alleles'], 1, [0,1], 'randomRestarts_A.svg')
	plotDataRestarts(noiseLevels, averagedMuErrors, groupedAboveStdMu, groupedBelowStdMu, ['Tumor fraction'], 3, [0,0.6], 'randomRestarts_Mu.svg')
	plotDataRestarts(noiseLevels, averagedTreeErrors, groupedAboveStdT, groupedBelowStdT, ['Trees'], 4, [-1,10], 'randomRestarts_T.svg')

def generateRandomRestartsFigure(dataFolder, noiseLevels):
	
	#F1. Read the data from all the simulation folders (The normal and permuted errors)
	[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors] = readData(dataFolder, noiseLevels, '')
	
	#F2. Average the errors per noise level
	[averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAmbiguityErrors, averagedSwapAbsentErrors, averagedSwapPresentErrors] = \
	averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAmbiguityErrors, groupedAncestryAbsentErrors, groupedAncestryPresentErrors)
	
	#Obtain the standard deviation above and below the mean for the averages
	[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
	
	[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
	[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
	[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)

	#F3. Plot the error per noise level in one figure
	
	plotFigureOneRandomRestarts(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors,
				  noiseLevels, groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu, groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT)
# Fig S6	
# dataFolder = 'Results/random_restarts'
# 
# noiseLevels = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
# generateRandomRestartsFigure(dataFolder, noiseLevels)
# exit()	

def plotAmbiguityScores(noiseLevels, ambiguities, ambiguityScores, ambiguityStds):
	
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	
	ax = plt.gca()
	
	legendLines = []
	labels = ['Ambiguities', 'Resolved ambiguities']
	p = ax.errorbar(noiseLevels, ambiguities, label='Ambiguities', color=colors[0], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevels, ambiguityScores, label='Resolved ambiguities', color=colors[1], linewidth=2, yerr=[ambiguityStds[1], ambiguityStds[0]])
	legendLines.append(p[0])
	
	ax.legend(legendLines, labels, loc=1, numpoints=1)
	
	ax.set_ylim(0,1)
	ax.set_xlabel(r'Sequencing noise ($\sigma$)')
 	ax.set_ylabel('Ratio')
	ax.set_xlim(-0.005,0.105)
	extraticks = [0.005, 0.01, 0.015, 0.025, 0.03]
	lim = ax.get_xlim()
	ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	ax.set_xticklabels(['', 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1, ''], rotation=90)
	ax.set_xlim(lim)
	
	plt.tight_layout()
	#plt.show()
	plt.savefig('ambiguities.svg')
	
			
#Figure 3G

#loop through the noise levels
# 
# noiseLevels = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1]
# #noiseLevels = [0, 0.005]
# ambiguities = []
# ambiguityScores = []
# 
# groupedAmbiguityScores = dict()
# 
# for noiseLevel in noiseLevels:
# 	
# 	currentSimulationFolder = simulationFolder + str(noiseLevel) + '/'
# 	print "noise level: ", noiseLevel
# 	[averageAmbiguities, averageAmbiguityScore, allAmbiguityScores] = computeCorrectAmbiguityScore(LAFAndCombinations, currentSimulationFolder)
# 	ambiguities.append(averageAmbiguities)
# 	ambiguityScores.append(averageAmbiguityScore)
# 	
# 	groupedAmbiguityScores[noiseLevel] = allAmbiguityScores
# 	
# #commpute the standard deviations
# ambiguityStds = obtainStandardDeviations(groupedAmbiguityScores, ambiguityScores)
# 
# plotAmbiguityScores(noiseLevels, ambiguities, ambiguityScores, ambiguityStds)

# exit()

#Plot the tree reconstruction errors when we use different measures. C vs. A, Euclidean, SNVs
#It is not fair to compare to the iterative scores, so we re-compute the score based on A using only one iteration. Otherwise the entire simulations will need to be re-run. 

#1. Read the A, C and SNV matrices from the simulation data folder

def readSimulationData(simulationFolderLocation, noiseLevels):
	
	allMethodTreeErrors = dict()
	allATreeErrors = dict()
	allCTreeErrors = dict()
	allEuclideanTreeErrors = dict()
	allSnvTreeErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = simulationFolderLocation + str(noiseLevel)
		
		methodTreeErrors = []
		aTreeErrors = []
		cTreeErrors = []
		euclideanTreeErrors = []
		snvTreeErrors = []
		
		for subdir, dirs, files in os.walk(simulationFolder):
			if subdir == simulationFolder: #we are not interested in the root folder
				continue
			for file in files:
				
				if re.match('treeError', file): #read the file and obtain the error
					methodTreeErrors += collectErrorsFromFile(file, subdir)
				if re.match('EstimatedA', file): #read the file and obtain the a matrix
					aMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
				if re.match('EstimatedC', file): #read the file and obtain the a matrix
					cMatrix = np.loadtxt(subdir + '/' + file, dtype=str)
				if re.match('LAFMeasurements_1', file): #read the file and obtain the error
					lafMatrix = np.loadtxt(subdir + '/' + file, dtype=float)
				if re.match('AFMeasurements_1', file): #read the file and obtain the error
					afMatrix = np.loadtxt(subdir + '/' + file, dtype=float)
				if re.match('RealTrees', file):
					stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
					tree = eval(stringDict)
					realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
				
				if re.match('simulationData.pkl', file):
					pkl_file = open(subdir + '/simulationData.pkl', 'rb')
					simulationData = pkl.load(pkl_file)
					pkl_file.close()
					
					#are the indices that we need in here?
					
					variantIndices = simulationData.samples[1].somaticVariantsInd
					chromosomes = simulationData.samples[1].measurements.chromosomes
					positions = simulationData.samples[1].measurements.starts
					segmentation = simulationData.samples[1].measurements.segmentation
						
					
					
				if re.match('SomVar', file):
					snvMatrix = np.loadtxt(subdir + '/' + file, dtype=int)
			
			snvTreeErrors.append(computeTreeErrorOtherMetrics.computeSNVTreeError(snvMatrix, cMatrix, lafMatrix, realTree, variantIndices, chromosomes, positions))
		
			#Generate a tree using the A matrix
			aTreeErrors.append(computeTreeErrorOtherMetrics.computeATreeError(aMatrix, lafMatrix, afMatrix, realTree, chromosomes, positions))
			
			#Generate a tree using the C matrix	
			cTreeErrors.append(computeTreeErrorOtherMetrics.computeCTreeError(cMatrix, realTree))
			
			#Generate a tree using the euclidean distance
			treeScore = computeTreeErrorOtherMetrics.computeEuclideanTreeError(lafMatrix, snvMatrix, realTree)
			euclideanTreeErrors.append(treeScore)
		
		#Compute errors based on SNVs?
		allATreeErrors[noiseLevel] = aTreeErrors
		allCTreeErrors[noiseLevel] = cTreeErrors
		allEuclideanTreeErrors[noiseLevel] = euclideanTreeErrors
		allSnvTreeErrors[noiseLevel] = snvTreeErrors
		allMethodTreeErrors[noiseLevel] = methodTreeErrors
	return [allATreeErrors, allCTreeErrors, allSnvTreeErrors, allEuclideanTreeErrors, allMethodTreeErrors]



def plotTreeErrorsDifferentMetrics(noiseLevels, averagedATreeErrors, averagedCTreeErrors, averagedSnvTreeErrors, averagedEuclideanTreeErrors, averagedMethodTreeErrors,
								   stdATreeErrors, stdCTreeErrors, stdSnvTreeErrors, stdEuclideanTreeErrors, stdMethodTreeErrors):
	
	
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	
	ax = plt.gca()
	
	#Define the offsets to remove overlap of error bars
	aOffset = 0.0003
	snvOffset = 0.0007
	euclOffset = 0.0011
	tcOffset = 0.00015
	
	noiseLevelsA = [i-aOffset for i in noiseLevels]
	noiseLevelsSnv = [i-snvOffset for i in noiseLevels]
	noiseLevelsEucl = [i-euclOffset for i in noiseLevels]
	noiseLevelsTc = [i-tcOffset for i in noiseLevels]
	
	#correct the standard deviations if these go below 0
	correctedBelowStd = []
	for std in range(0, len(stdMethodTreeErrors[1])):
		newStd = stdMethodTreeErrors[1][std]
		if (averagedMethodTreeErrors[std]-newStd) < 0:
			newStd = abs(0-averagedMethodTreeErrors[std])
		correctedBelowStd.append(newStd)
	
	
	legendLines = []
	labels = ['Alleles', 'Copy numbers', 'Somatic SNV', 'Euclidean (LAF)', 'TargetClone']
	p = ax.errorbar(noiseLevelsA, averagedATreeErrors, yerr=[stdATreeErrors[1], stdATreeErrors[0]], label='A', color=colors[1], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevels, averagedCTreeErrors, yerr=[stdCTreeErrors[1], stdCTreeErrors[0]], label='C', color=colors[0], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevelsSnv, averagedSnvTreeErrors, yerr=[stdSnvTreeErrors[1], stdSnvTreeErrors[0]], label='SNV', color=colors[3], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevelsEucl, averagedEuclideanTreeErrors, yerr=[stdEuclideanTreeErrors[1], stdEuclideanTreeErrors[0]], label='Euclidean', color=colors[4], linewidth=2)
	legendLines.append(p[0])
	p = ax.errorbar(noiseLevelsTc, averagedMethodTreeErrors, yerr=[correctedBelowStd, stdMethodTreeErrors[0]], label='TargetClone', color=colors[5], linewidth=2)
	legendLines.append(p[0])
	
	ax.legend(legendLines, labels, loc=2, numpoints=1)
	
	ax.set_ylim(-0.5,8)
	ax.set_xlabel(r'Sequencing noise ($\sigma$)')
 	ax.set_ylabel('Error (trees)')
	#ax.set_xlim(-0.005,0.105)
	ax.set_xlim(0,0.035)
	#extraticks = [0.005, 0.01, 0.015, 0.025, 0.03]
	lim = ax.get_xlim()
	#ax.set_xticks(sorted(list(ax.get_xticks()) + extraticks))
	#ax.set_xticklabels(['', 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08, 0.1, ''], rotation=90)
	ax.set_xticklabels(['', 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, ''], rotation=90)
	ax.set_xlim(lim)
	
	plt.tight_layout()
	#plt.show()
	plt.savefig('CvsAvsSNVvsEuclidean.svg')
	

#Figure S7
# 
# print "parsing simulation data"
# #Obtain all errors
# simulationFolder = 'Results/generic_noise'
# noiseLevels = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
# [allATreeErrors, allCTreeErrors, allSnvTreeErrors, allEuclideanTreeErrors, allMethodTreeErrors] = readSimulationData(simulationFolder, noiseLevels)
# #Compute the average errors and standard deviations
# print "computing averages: "
# averagedATreeErrors = averageData(allATreeErrors, 'T')
# averagedCTreeErrors = averageData(allCTreeErrors, 'T')
# averagedSnvTreeErrors = averageData(allSnvTreeErrors, 'T')
# averagedEuclideanTreeErrors = averageData(allEuclideanTreeErrors, 'T')
# averagedMethodTreeErrors = averageData(allMethodTreeErrors, 'T')
# print "computing standard deviations: "
# stdATreeErrors = obtainStandardDeviations(allATreeErrors, averagedATreeErrors)
# print "c:"
# stdCTreeErrors = obtainStandardDeviations(allCTreeErrors, averagedCTreeErrors)
# print "snv:"
# stdSnvTreeErrors = obtainStandardDeviations(allSnvTreeErrors, averagedSnvTreeErrors)
# print "euclidean:"
# stdEuclideanTreeErrors = obtainStandardDeviations(allEuclideanTreeErrors, averagedEuclideanTreeErrors)
# print "tree:"
# stdMethodTreeErrors = obtainStandardDeviations(allMethodTreeErrors, averagedMethodTreeErrors)
# 
# allATreeErrors = None
# 
# 
# print "plotting: "
# #Make a plot of the tree errors
# plotTreeErrorsDifferentMetrics(noiseLevels, averagedATreeErrors, averagedCTreeErrors, averagedSnvTreeErrors, averagedEuclideanTreeErrors, averagedMethodTreeErrors, stdATreeErrors, stdCTreeErrors, stdSnvTreeErrors,
# 							   stdEuclideanTreeErrors, stdMethodTreeErrors)
# 
# exit()






#This makes bins for C, A, mu and T separately
#THese bins can then be provided to the plotting function
def binData(sortedLeftOut, sortedErrors, sortedAmbiguityCorrectedErrors):
	binnedLeftOut = []
	binnedError = []
	binnedStdAbove = []
	binnedStdBelow = []
	binnedErrorP = []
	binnedAmbiguityErrors = []
	
	previousLeftOut = 1
	currentErrors = []
	currentAmbiguityErrors = []
	binnedLeftOut.append(sortedLeftOut[0]) #add the first one
	for leftOut in range(0, len(sortedLeftOut)):
		
		#if we find a new left out value, take the mean of all current values and start anew
		if sortedLeftOut[leftOut] != previousLeftOut:
			meanError = sum(currentErrors) / float(len(currentErrors))
			meanAmbiguityError = sum(currentAmbiguityErrors) / float(len(currentAmbiguityErrors))
			binnedError.append(meanError)
			binnedAmbiguityErrors.append(meanAmbiguityError)
			
			#Check how many of the values are above/below the mean for the error bars
			above = []
			below = []
			
			for value in range(0, len(currentErrors)):
				if currentErrors[value] > meanError:
					above.append(currentErrors[value])
				else:
					below.append(currentErrors[value])
			
			binnedStdAbove.append(above)
			binnedStdBelow.append(below)
			
			binnedLeftOut.append(sortedLeftOut[leftOut])
			currentErrors = []
			currentAmbiguityErrors = []
		else:
			currentErrors.append(sortedErrors[leftOut])
			currentAmbiguityErrors.append(sortedAmbiguityCorrectedErrors[leftOut])
		previousLeftOut = sortedLeftOut[leftOut]
		
		
	#Also add the last values
	meanError = sum(currentErrors) / float(len(currentErrors))
	meanAmbiguityError = sum(currentAmbiguityErrors) / float(len(currentAmbiguityErrors))
	binnedError.append(meanError)
	binnedAmbiguityErrors.append(meanAmbiguityError)
	
	#Check how many of the values are above/below the mean for the error bars
	above = []
	below = []
	
	for value in range(0, len(currentErrors)):
		if currentErrors[value] > meanError:
			above.append(currentErrors[value])
		else:
			below.append(currentErrors[value])
	
	binnedStdAbove.append(np.std(above))
	binnedStdBelow.append(np.std(below))
	
	return [binnedLeftOut, binnedError, binnedStdAbove, binnedStdBelow]
# 
# [binnedLeftOut, binnedC, binnedStdAboveC, binnedStdBelowC] = binData(sortedLeftOut, sortedCErrors, sortedAmbiguityErrors) #exclude the ambiguities for now, not interesting
# [binnedLeftOut, binnedA, binnedStdAboveA, binnedStdBelowA] = binData(sortedLeftOut, sortedAErrors, sortedAmbiguityErrors)
# [binnedLeftOut, binnedMu, binnedStdAboveMu, binnedStdBelowMu] = binData(sortedLeftOut, sortedMuErrors, sortedAmbiguityErrors)
# [binnedLeftOut, binnedT, binnedStdAboveT, binnedStdBelowT] = binData(sortedLeftOut, sortedTErrors, sortedAmbiguityErrors)
# 










def sortData(simulationFolder):
	
	orderedMu = []
	orderedCErrors = []
	orderedAErrors = []
	orderedMuErrors = []
	orderedTErrors = []
	
	ambiguityErrors = []
	ambiguityCorrectedErrors = []
	
	#Repeat for the permutations

	orderedPCErrors = []
	orderedPAErrors = []
	orderedPMuErrors = []
	orderedPTErrors = []
	
	ambiguityPErrors = []
	ambiguityPCorrectedErrors = []
	
	#This needs an ordering
	for subdir, dirs, files in os.walk(simulationFolder):
		#print subdir
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		currentMu = 0
		cError = 0
		aError = 0
		muError = 0
		tError = 0
		ambiguityError = 0
		ambiguityCorrectedError = 0
		
		pCError = 0
		pAError = 0
		pMuError = 0
		pTError = 0
		pAmbiguityError = 0
		pAmbiguityCorrectedError = 0
		
		#print subdir
		
		for file in files:
			#print file
			if re.match('RealMu', file): #read the file and obtain the error
				
				currentMu = float(collectErrorsFromFile(file, subdir)[1])
				
			if re.match('cError', file):
				cError = float(collectErrorsFromFile(file, subdir)[0])

			if re.match('aError', file):
				aError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('muError', file):
				muError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('treeError', file):
				tError = float(collectErrorsFromFile(file, subdir)[0])
			
			if re.match('ambiguityError', file):
				ambiguityError = float(collectErrorsFromFile(file, subdir)[0])
			if re.match('ambiguityCorrectedError', file):
				ambiguityCorrectedError = float(collectErrorsFromFile(file, subdir)[0])
				
			
			if re.match('pCError', file):
				pCError = float(collectErrorsFromFile(file, subdir)[0])

			if re.match('pAError', file):
				pAError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('pMuError', file):
				pMuError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('pTreeError', file):
				pTError = float(collectErrorsFromFile(file, subdir)[0])
			
			if re.match('pAmbiguityError', file):
				pAmbiguityError = float(collectErrorsFromFile(file, subdir)[0])
			if re.match('pAmbiguityCorrectedError', file):
				pAmbiguityCorrectedError = float(collectErrorsFromFile(file, subdir)[0])	
		
		
		orderedMu.append(currentMu)
		orderedCErrors.append(cError)
		orderedAErrors.append(aError)
		orderedMuErrors.append(muError)
		orderedTErrors.append(tError)
		ambiguityErrors.append(ambiguityError)
		ambiguityCorrectedErrors.append(ambiguityCorrectedError)
		
		orderedPCErrors.append(pCError)
		orderedPAErrors.append(pAError)
		orderedPMuErrors.append(pMuError)
		orderedPTErrors.append(pTError)
		ambiguityPErrors.append(pAmbiguityError)
		ambiguityPCorrectedErrors.append(pAmbiguityCorrectedError)
		
	#Then sort by mu
		
	res = [i[0] for i in sorted(enumerate(orderedMu), key=lambda x:x[1])]
	
	sortedMu = []
	sortedCErrors = []
	sortedAErrors = []
	sortedMuErrors = []
	sortedTErrors = []
	sortedAmbiguityErrors = []
	sortedAmbiguityCorrectedErrors = []
	
	sortedPCErrors = []
	sortedPAErrors = []
	sortedPMuErrors = []
	sortedPTErrors = []
	sortedPAmbiguityErrors = []
	sortedPAmbiguityCorrectedErrors = []
	
	for ind in res:
		sortedMu.append(orderedMu[ind])
		sortedCErrors.append(orderedCErrors[ind])
		sortedAErrors.append(orderedAErrors[ind])
		sortedMuErrors.append(orderedMuErrors[ind])
		sortedTErrors.append(orderedTErrors[ind])
		sortedAmbiguityErrors.append(ambiguityErrors[ind])
		sortedAmbiguityCorrectedErrors.append(ambiguityCorrectedErrors[ind])
		
		sortedPCErrors.append(orderedPCErrors[ind])
		sortedPAErrors.append(orderedPAErrors[ind])
		sortedPMuErrors.append(orderedPMuErrors[ind])
		sortedPTErrors.append(orderedPTErrors[ind])
		sortedPAmbiguityErrors.append(ambiguityPErrors[ind])
		sortedPAmbiguityCorrectedErrors.append(ambiguityPCorrectedErrors[ind])
	
	
	return [sortedMu, sortedCErrors, sortedAErrors, sortedMuErrors, sortedTErrors, sortedAmbiguityErrors, sortedAmbiguityCorrectedErrors, sortedPCErrors, sortedPAErrors, sortedPMuErrors, sortedPTErrors, sortedPAmbiguityErrors, sortedPAmbiguityCorrectedErrors]
	

#Call the function for all noise levels
#Have a function that bins by mu

def binValues(orderedMu, sortedErrors):
	binnedMu = []
	binnedError = []
	binnedStdAbove = []
	binnedStdBelow = []
	binnedAmbiguityErrors = []
	
	binSize = 10
	boxCount = len(orderedMu) / binSize
	offset = 0
	for boxInd in range(0, boxCount):
		
		errorValues = [sortedErrors[i] for i in range(offset,offset+binSize)]
		binnedMu.append(boxInd)
		meanError = sum(errorValues) / float(len(errorValues))
		binnedError.append(meanError)
	
		
		#Compute all values above and below the mean
		above = []
		below = []
		
		
		for value in range(0, len(errorValues)):
			if errorValues[value] > meanError:
				above.append(errorValues[value])
			else:
				below.append(errorValues[value])
			

		binnedStdAbove.append(np.std(errorValues))
		binnedStdBelow.append(np.std(errorValues))
		
		#Make sure that the sd is also reported per bin
		
		
		offset += binSize

	return [binnedMu, binnedError, binnedStdAbove, binnedStdBelow]
# 
# 
# #Bin all values that we are interested in
# [binnedMu, binnedCErrors0, binnedStdAboveC0, binnedStdBelowC0] = binValues(sortedMu0, sortedCErrors0)
# [binnedMu, binnedCErrors2, binnedStdAboveC2, binnedStdBelowC2]  = binValues(sortedMu2, sortedCErrors2)
# [binnedMu, binnedCErrors4, binnedStdAboveC4, binnedStdBelowC4]  = binValues(sortedMu4, sortedCErrors4)
# [binnedMu, binnedCErrors6, binnedStdAboveC6, binnedStdBelowC6]  = binValues(sortedMu6, sortedCErrors6)
# [binnedMu, binnedCErrors8, binnedStdAboveC8, binnedStdBelowC8]  = binValues(sortedMu8, sortedCErrors8)
# [binnedMu, binnedCErrors10, binnedStdAboveC10, binnedStdBelowC10]  = binValues(sortedMu10, sortedCErrors10)
# 
# [binnedMu, binnedAErrors0, binnedStdAboveA0, binnedStdBelowA0] = binValues(sortedMu0, sortedAErrors0)
# [binnedMu, binnedAErrors2, binnedStdAboveA2, binnedStdBelowA2]  = binValues(sortedMu2, sortedAErrors2)
# [binnedMu, binnedAErrors4, binnedStdAboveA4, binnedStdBelowA4]  = binValues(sortedMu4, sortedAErrors4)
# [binnedMu, binnedAErrors6, binnedStdAboveA6, binnedStdBelowA6]  = binValues(sortedMu6, sortedAErrors6)
# [binnedMu, binnedAErrors8, binnedStdAboveA8, binnedStdBelowA8]  = binValues(sortedMu8, sortedAErrors8)
# [binnedMu, binnedAErrors10, binnedStdAboveA10, binnedStdBelowA10]  = binValues(sortedMu10, sortedAErrors10)
# 
# [binnedMu, binnedMuErrors0, binnedStdAboveMu0, binnedStdBelowMu0] = binValues(sortedMu0, sortedMuErrors0)
# [binnedMu, binnedMuErrors2, binnedStdAboveMu2, binnedStdBelowMu2]  = binValues(sortedMu2, sortedMuErrors2)
# [binnedMu, binnedMuErrors4, binnedStdAboveMu4, binnedStdBelowMu4]  = binValues(sortedMu4, sortedMuErrors4)
# [binnedMu, binnedMuErrors6, binnedStdAboveMu6, binnedStdBelowMu6]  = binValues(sortedMu6, sortedMuErrors6)
# [binnedMu, binnedMuErrors8, binnedStdAboveMu8, binnedStdBelowMu8]  = binValues(sortedMu8, sortedMuErrors8)
# [binnedMu, binnedMuErrors10, binnedStdAboveMu10, binnedStdBelowMu10]  = binValues(sortedMu10, sortedMuErrors10)
# 
# [binnedMu, binnedTErrors0, binnedStdAboveT0, binnedStdBelowT0] = binValues(sortedMu0, sortedTErrors0)
# [binnedMu, binnedTErrors2, binnedStdAboveT2, binnedStdBelowT2]  = binValues(sortedMu2, sortedTErrors2)
# [binnedMu, binnedTErrors4, binnedStdAboveT4, binnedStdBelowT4]  = binValues(sortedMu4, sortedTErrors4)
# [binnedMu, binnedTErrors6, binnedStdAboveT6, binnedStdBelowT6]  = binValues(sortedMu6, sortedTErrors6)
# [binnedMu, binnedTErrors8, binnedStdAboveT8, binnedStdBelowT8]  = binValues(sortedMu8, sortedTErrors8)
# [binnedMu, binnedTErrors10, binnedStdAboveT10, binnedStdBelowT10]  = binValues(sortedMu10, sortedTErrors10)
# 
# #Do the same binning but then for the completely random case
# 
# 
# 
# 
# #Group the values
# binnedCErrors = [binnedCErrors0, binnedCErrors2, binnedCErrors4, binnedCErrors6, binnedCErrors8, binnedCErrors10]
# binnedCStdAbove = [binnedStdAboveC0, binnedStdAboveC2, binnedStdAboveC4, binnedStdAboveC6, binnedStdAboveC8, binnedStdAboveC10]
# binnedCStdBelow = [binnedStdBelowC0, binnedStdBelowC2, binnedStdBelowC4, binnedStdBelowC6, binnedStdBelowC8, binnedStdBelowC10]
# 
# binnedAErrors = [binnedAErrors0, binnedAErrors2, binnedAErrors4, binnedAErrors6, binnedAErrors8, binnedAErrors10]
# binnedAStdAbove = [binnedStdAboveA0, binnedStdAboveA2, binnedStdAboveA4, binnedStdAboveA6, binnedStdAboveA8, binnedStdAboveA10]
# binnedAStdBelow = [binnedStdBelowA0, binnedStdBelowA2, binnedStdBelowA4, binnedStdBelowA6, binnedStdBelowA8, binnedStdBelowA10]
# 
# binnedMuErrors = [binnedMuErrors0, binnedMuErrors2, binnedMuErrors4, binnedMuErrors6, binnedMuErrors8, binnedMuErrors10]
# binnedMuStdAbove = [binnedStdAboveMu0, binnedStdAboveMu2, binnedStdAboveMu4, binnedStdAboveMu6, binnedStdAboveMu8, binnedStdAboveMu10]
# binnedMuStdBelow = [binnedStdBelowMu0, binnedStdBelowMu2, binnedStdBelowMu4, binnedStdBelowMu6, binnedStdBelowMu8, binnedStdBelowMu10]
# 
# binnedTErrors = [binnedTErrors0, binnedTErrors2, binnedTErrors4, binnedTErrors6, binnedTErrors8, binnedTErrors10]
# binnedTStdAbove = [binnedStdAboveT0, binnedStdAboveT2, binnedStdAboveT4, binnedStdAboveT6, binnedStdAboveT8, binnedStdAboveT10]
# binnedTStdBelow = [binnedStdBelowT0, binnedStdBelowT2, binnedStdBelowT4, binnedStdBelowT6, binnedStdBelowT8, binnedStdBelowT10]



# # 
simulationFolder = 'Results/generic_noise0.02/'
[sortedMu2, sortedCErrors2, sortedAErrors2, sortedMuErrors2, sortedTErrors2, sortedAmbiguityScores2, sortedAmbiguityCorrectedErrors2,
 sortedPCErrors2, sortedPAErrors2, sortedPMuErrors2, sortedPTErrors2, sortedPAmbiguityErrors2, sortedPAmbiguityCorrectedErrors2] = sortData(simulationFolder)

[binnedMu, binnedCErrors2, binnedStdAboveC2, binnedStdBelowC2]  = binValues(sortedMu2, sortedCErrors2)

[binnedMu, binnedAErrors2, binnedStdAboveA2, binnedStdBelowA2]  = binValues(sortedMu2, sortedAErrors2)

[binnedMu, binnedMuErrors2, binnedStdAboveMu2, binnedStdBelowMu2]  = binValues(sortedMu2, sortedMuErrors2)

[binnedMu, binnedTErrors2, binnedStdAboveT2, binnedStdBelowT2]  = binValues(sortedMu2, sortedTErrors2)


#Group the values
binnedCErrors = [binnedCErrors2]
binnedCStd = [binnedStdAboveC2, binnedStdBelowC2]

binnedAErrors = [binnedAErrors2]
binnedAStd = [binnedStdAboveA2, binnedStdBelowA2]

binnedMuErrors = [binnedMuErrors2]
binnedMuStd = [binnedStdAboveMu2, binnedStdBelowMu2]

binnedTErrors = [binnedTErrors2]
binnedTStd = [binnedStdAboveT2, binnedStdBelowT2]

#Plot the noise level values in the same figure for C, A, mu and T

def plotNoise2TumorFractionsSummary(mu, binnedC, binnedA, binnedMu, binnedT, binnedCStd, binnedAStd, binnedMuStd, binnedTStd):
	
	plt.figure()
	colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	#Group the errors in one figure for C, A, Mu and T

	#make the plot
	aOffset = 0.07
	muOffset = 0.11
	treeOffset = 0.17
	
	noiseLevelsA = [i-aOffset for i in mu]
	noiseLevelsMu = [i-muOffset for i in mu]
	noiseLevelsTree = [i-treeOffset for i in mu]
	print mu

	ax = axes()
	ax2 = ax.twinx()
	
	labels = []
	legendLines = []
	legendLines2 = []
	labels2 = []
	
	print "plotting data: "
	print mu
	print binnedC[0]
	print binnedCStd[1]
	print binnedCStd[0]
	
	
	stdBelow = []
	for error in range(0, len(binnedTStd[1])):
		newStd = binnedTStd[1][error]
		if binnedT[0][error]-newStd < 0:
			newStd = abs(0-binnedT[0][error])
		stdBelow.append(newStd)
		
	
	p = ax.errorbar(mu, binnedC[0], yerr=[binnedCStd[1], binnedCStd[0]], label="C error", color = colors[0], linewidth=2, antialiased=True)
	legendLines.append(p[0])
	labels.append('Copy numbers')
	p = ax.errorbar(noiseLevelsA, binnedA[0], yerr=[binnedAStd[1], binnedAStd[0]], label="A error", color = colors[1], linewidth=2, antialiased=True)
	legendLines.append(p[0])
	labels.append('Alleles')
	p = ax.errorbar(noiseLevelsMu, binnedMu[0], yerr=[binnedMuStd[1], binnedMuStd[0]], label= "Mu error", color = colors[3], linewidth=2, antialiased=True)
	legendLines.append(p[0])
	labels.append('Tumor fraction')
	p = ax2.errorbar(noiseLevelsTree, binnedT[0], yerr=[stdBelow, binnedTStd[0]], label="T error", color = colors[4], linewidth=2, antialiased=True)
	legendLines2.append(p[0])
	labels2.append('Trees')

	
	ax.set_ylim(0,0.8)
	ax2.set_ylim(-0.5,12)
	
	ax.set_xlim(-1,10)
	
	ax.legend(legendLines, labels, loc=2)
	ax2.legend(legendLines2, labels2)
	ax.set_xlabel('mu (% tumor fraction)')
 	ax.set_ylabel('Error')
	ax2.set_ylabel('Error')
	#ax.set_xticklabels(['', 1, 2, 3, 4, 5, 6, ''])
	ax.set_xticklabels(['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100'], rotation=90)
	ax.set_xticks(range(0,10))
	ax.set_xlabel('Mu (% tumor fraction)')
	ax.set_ylabel('Error')
	#splt.show()
	plt.tight_layout()
	plt.savefig('tumorFractionSummary.svg')

#Figure 3F
# 	
# plotNoise2TumorFractionsSummary(binnedMu, binnedCErrors, binnedAErrors, binnedMuErrors, binnedTErrors, binnedCStd, binnedAStd, binnedMuStd, binnedTStd)
# 
# exit()


#Make scatterplots showing the difference between the normal simulations and the horizontal permutations




def binPermutedData(permutedData):
	print permutedData
	noiseLevel = permutedData.keys()[0]
	#bin every 10 positions
	binnedPermutedData = []
	binSize = 10
	currentValues = []
	offset = 0
	for i in range(0,1001):
		print "i is : ", i+offset
		if i < (offset+binSize)-1:
			currentValues.append(permutedData[noiseLevel][i])
			
		else:
			offset += binSize
			currentAverage = sum(currentValues)/binSize
			binnedPermutedData.append(currentAverage)
			currentValues = []
			i+=offset
			
	currentAverage = sum(currentValues)/binSize
	binnedPermutedData.append(currentAverage)
	
	return binnedPermutedData
from scipy import stats
def density_estimation(m1, m2):
    X, Y = np.mgrid[min(m1)-1:max(m1)+1:100j, min(m2)-1:max(m2)+1:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z	

def plotNormalPermutedComparison(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors):
	
	from scipy.stats import gaussian_kde
	
	#plt.figure()
	ax = axes()
	
	x = groupedCErrors
	y = groupedPCErrors
	
	# Calculate the point density
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	
	plt.scatter(groupedCErrors, groupedPCErrors, c=z, edgecolor='')
	X, Y, Z = density_estimation(groupedCErrors, groupedPCErrors)
	plt.contour(X, Y, Z)
	ax.set_xlim([min(groupedCErrors)-0.025, max(groupedCErrors)+0.025])                                                                           
	ax.set_ylim([min(groupedPCErrors)-0.075, max(groupedPCErrors)+0.025])    
	ax.set_xlabel('Simulation (error)')
	ax.set_ylabel('Horizontal permutation (error)')
	plt.plot( groupedCErrors,groupedCErrors, antialiased=True ) # identity line
	plt.savefig('C_horizontalPermutations_rerun2.svg')
	#plt.show()
	#exit()
	
	plt.figure()
	ax = axes()
	
	x = groupedAErrors
	y = groupedPAErrors
	
	# Calculate the point density
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	plt.scatter(groupedAErrors, groupedPAErrors, c=z, edgecolor='')
	X, Y, Z = density_estimation(groupedAErrors, groupedPAErrors)
	plt.contour(X, Y, Z)
	ax.set_xlim([min(groupedAErrors)-0.025, max(groupedAErrors)+0.025])                                                                           
	ax.set_ylim([min(groupedPAErrors)-0.075, max(groupedPAErrors)+0.025])    
	ax.set_xlabel('Simulation (error)')
	ax.set_ylabel('Horizontal permutation (error)')
	plt.plot( groupedAErrors,groupedAErrors ) # identity line
	plt.savefig('A_horizontalPermutations_rerun2.svg')
	#plt.show()
	
	plt.figure()
	ax = axes()
	
	x = groupedMuErrors
	y = groupedPMuErrors
	
	# Calculate the point density
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	plt.scatter(groupedMuErrors, groupedPMuErrors, c=z, edgecolor='')
	
	#Add contour lines
	X, Y, Z = density_estimation(groupedMuErrors, groupedPMuErrors)
	plt.contour(X, Y, Z)
	ax.set_xlim([min(groupedMuErrors)-0.01, max(groupedMuErrors)+0.010])                                                                           
	ax.set_ylim([min(groupedPMuErrors)-0.00, max(groupedPMuErrors)+0.03])    
	ax.set_xlabel('Simulation (error)')
	ax.set_ylabel('Horizontal permutation (error)')
	plt.plot( groupedMuErrors,groupedMuErrors ) # identity line
	plt.savefig('Mu_horizontalPermutations_rerun2.svg')
	#plt.show()
	
	#Calculate the mean and standard deviation for both axes
	meanX = sum(groupedTreeErrors)/float(len(groupedTreeErrors))
	meanY = sum(groupedPTreeErrors)/float(len(groupedPTreeErrors))
	# 
	# pointsAbove = []
	# pointsBelow = []
	# 
	# for error in groupedTreeErrors:
	# 	if error >= meanX:
	# 		pointsAbove.append(error)
	# 	else:
	# 		pointsBelow.append(error)
	# 		
	# stdAboveX = np.std(pointsAbove)
	# stdBelowX = np.std(pointsBelow)
	# from math import sqrt
	# 
	# differences = [x - meanX for x in pointsAbove]
	# sq_differences = [d ** 2 for d in differences]
	# ssd = sum(sq_differences)
	# variance = ssd / (len(pointsAbove) - 1)
	# stdAboveX = sqrt(variance)
	# 
	# differences = [x - meanX for x in pointsBelow]
	# sq_differences = [d ** 2 for d in differences]
	# ssd = sum(sq_differences)
	# variance = ssd / (len(pointsBelow) - 1)
	# stdBelowX = sqrt(variance)
	# 
	# 
	# 
	# pointsAbove = []
	# pointsBelow = []
	# 
	# for error in groupedPTreeErrors:
	# 	if error >= meanY:
	# 		pointsAbove.append(error)
	# 	else:
	# 		pointsBelow.append(error)
	# 		
	# differences = [x - meanY for x in pointsAbove]
	# sq_differences = [d ** 2 for d in differences]
	# ssd = sum(sq_differences)
	# variance = ssd / (len(pointsAbove) - 1)
	# stdAboveY = sqrt(variance)
	# 
	# differences = [x - meanY for x in pointsBelow]
	# sq_differences = [d ** 2 for d in differences]
	# ssd = sum(sq_differences)
	# variance = ssd / (len(pointsBelow) - 1)
	# stdBelowY = sqrt(variance)
	# 
	x = groupedTreeErrors
	y = groupedPTreeErrors
	
	
	
	#Add jitter
	#For every point, we add a random offset (either + or -)
	from random import randint
	
	jitterX = []
	
	for error in x:
	
		randOffset = randint(0,25) /float(100)
		randDirectionality = randint(0,1)
		if randDirectionality > 0:
			randOffset = -randOffset
		jitterX.append(error+randOffset)
	
	jitterY = []
	
	for error in y:
	
		randOffset = randint(0,25) /float(100)
		randDirectionality = randint(0,1)
		if randDirectionality > 0:
			randOffset = -randOffset
		jitterY.append(error+randOffset)
	
	X, Y, Z = density_estimation(jitterX, jitterY)

	fig, ax = plt.subplots()                   
	xy = np.vstack([jitterX,jitterY])
	z = gaussian_kde(xy)(xy)
	
	plt.scatter(jitterX, jitterY, c=z, edgecolor='')
	# Show density 
	#ax.imshow(np.rot90(Z),                                                    
	#		  extent=[min(jitterX), max(jitterX), min(jitterY), max(jitterY)])
	
	# Add contour lines
	plt.contour(X, Y, Z)
	
	#error bars , , , yerr=[stdBelowY, stdAboveY]
	#ax.errorbar(meanX, meanY, antialiased=True)
	#ax.errorbar(meanX, meanY, yerr=[[stdBelowY], [stdAboveY]], xerr=[[stdBelowX], [stdAboveX]], linewidth=1, color='k', antialiased=True)
	#ax.plot(jitterX, jitterY, 'k.', markersize=2)    
	plt.plot( groupedTreeErrors,groupedTreeErrors, color='b') # identity line
	ax.set_xlim([min(jitterX)-1, max(jitterX)+1])                                                                           
	ax.set_ylim([min(jitterY)-1, max(jitterY)+1])                                                                           
	#plt.show()
	
	ax.set_xlabel('Simulation (error)')
	ax.set_ylabel('Horizontal permutation (error)')
	
	plt.savefig('T_horizontalPermutations_rerun2.svg')
	#plt.show()
	exit()
	
#Figure S8

#In the read data with permutations I added Mu75 after the path for now, quick and dirty
dataFolder = 'Results/generic_noise'
noiseLevels = ['0.02']
[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedPCErrors, groupedPAErrors, groupedPMuErrors, groupedPTreeErrors, a, b] = readDataIncludingPermutations(dataFolder, noiseLevels, '')

# #average the permuted data in bins
# #averagedPCErrors = binPermutedData(groupedPCErrors)
# #averagedPAErrors = binPermutedData(groupedPAErrors)
# #averagedPMuErrors = binPermutedData(groupedPMuErrors)
# #averagedPTreeErrors = binPermutedData(groupedPTreeErrors)
# 
# averagedPCErrors = groupedPCErrors
# averagedPAErrors = groupedPAErrors
# averagedPMuErrors = groupedPMuErrors
# averagedPTreeErrors = groupedPTreeErrors
# 
# 
# #plotNormalPermutedComparison(groupedCErrors[0.02], groupedAErrors[0.02], groupedMuErrors[0.02], groupedTreeErrors[0.02], averagedPCErrors, averagedPAErrors, averagedPMuErrors, averagedPTreeErrors)
# plotNormalPermutedComparison(groupedCErrors[0.02], groupedAErrors[0.02], groupedMuErrors[0.02], groupedTreeErrors[0.02], averagedPCErrors[0.02], averagedPAErrors[0.02], averagedPMuErrors[0.02], averagedPTreeErrors[0.02])
# 
# exit()

#Figure S5

#Also make a plot showing the correlation between A and T (at one noise level as well)

def plotATCorrelation(groupedAErrors, groupedTErrors):
	
	plt.figure()
	ax = axes()
	plt.scatter(groupedAErrors, groupedTErrors)
	ax.set_xlabel('Alleles (error)')
	ax.set_ylabel('Trees (error)')
	plt.savefig('A_TCorrelation.svg')

# 
# plotATCorrelation(groupedAErrors['0.02'], groupedTreeErrors['0.02'])
# exit()

#Figure S9

def sortData(simulationFolder):
	
	orderedMu = []
	orderedCErrors = []
	orderedAErrors = []
	orderedMuErrors = []
	orderedTErrors = []
	
	ambiguityErrors = []
	ambiguityCorrectedErrors = []
	
	#Repeat for the permutations

	orderedPCErrors = []
	orderedPAErrors = []
	orderedPMuErrors = []
	orderedPTErrors = []
	
	ambiguityPErrors = []
	ambiguityPCorrectedErrors = []
	
	#This needs an ordering
	for subdir, dirs, files in os.walk(simulationFolder):
		#print subdir
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		currentMu = 0
		cError = 0
		aError = 0
		muError = 0
		tError = 0
		ambiguityError = 0
		ambiguityCorrectedError = 0
		
		pCError = 0
		pAError = 0
		pMuError = 0
		pTError = 0
		pAmbiguityError = 0
		pAmbiguityCorrectedError = 0
		
		#print subdir
		
		for file in files:
			#print file
			if re.match('RealMu', file): #read the file and obtain the error
				
				currentMu = float(collectErrorsFromFile(file, subdir)[1])
				
			if re.match('cError', file):
				cError = float(collectErrorsFromFile(file, subdir)[0])

			if re.match('aError', file):
				aError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('muError', file):
				muError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('treeError', file):
				tError = float(collectErrorsFromFile(file, subdir)[0])
			
			if re.match('ambiguityError', file):
				ambiguityError = float(collectErrorsFromFile(file, subdir)[0])
			if re.match('ambiguityCorrectedError', file):
				ambiguityCorrectedError = float(collectErrorsFromFile(file, subdir)[0])
				
			
			if re.match('pCError', file):
				pCError = float(collectErrorsFromFile(file, subdir)[0])

			if re.match('pAError', file):
				pAError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('pMuError', file):
				pMuError = float(collectErrorsFromFile(file, subdir)[0])
				
			if re.match('pTreeError', file):
				pTError = float(collectErrorsFromFile(file, subdir)[0])
			
			if re.match('pAmbiguityError', file):
				pAmbiguityError = float(collectErrorsFromFile(file, subdir)[0])
			if re.match('pAmbiguityCorrectedError', file):
				pAmbiguityCorrectedError = float(collectErrorsFromFile(file, subdir)[0])	
		
		
		orderedMu.append(currentMu)
		orderedCErrors.append(cError)
		orderedAErrors.append(aError)
		orderedMuErrors.append(muError)
		orderedTErrors.append(tError)
		ambiguityErrors.append(ambiguityError)
		ambiguityCorrectedErrors.append(ambiguityCorrectedError)
		
		orderedPCErrors.append(pCError)
		orderedPAErrors.append(pAError)
		orderedPMuErrors.append(pMuError)
		orderedPTErrors.append(pTError)
		ambiguityPErrors.append(pAmbiguityError)
		ambiguityPCorrectedErrors.append(pAmbiguityCorrectedError)
		
	#Then sort by mu
		
	res = [i[0] for i in sorted(enumerate(orderedMu), key=lambda x:x[1])]
	
	sortedMu = []
	sortedCErrors = []
	sortedAErrors = []
	sortedMuErrors = []
	sortedTErrors = []
	sortedAmbiguityErrors = []
	sortedAmbiguityCorrectedErrors = []
	
	sortedPCErrors = []
	sortedPAErrors = []
	sortedPMuErrors = []
	sortedPTErrors = []
	sortedPAmbiguityErrors = []
	sortedPAmbiguityCorrectedErrors = []
	
	for ind in res:
		sortedMu.append(orderedMu[ind])
		sortedCErrors.append(orderedCErrors[ind])
		sortedAErrors.append(orderedAErrors[ind])
		sortedMuErrors.append(orderedMuErrors[ind])
		sortedTErrors.append(orderedTErrors[ind])
		sortedAmbiguityErrors.append(ambiguityErrors[ind])
		sortedAmbiguityCorrectedErrors.append(ambiguityCorrectedErrors[ind])
		
		sortedPCErrors.append(orderedPCErrors[ind])
		sortedPAErrors.append(orderedPAErrors[ind])
		sortedPMuErrors.append(orderedPMuErrors[ind])
		sortedPTErrors.append(orderedPTErrors[ind])
		sortedPAmbiguityErrors.append(ambiguityPErrors[ind])
		sortedPAmbiguityCorrectedErrors.append(ambiguityPCorrectedErrors[ind])
	
	
	return [sortedMu, sortedCErrors, sortedAErrors, sortedMuErrors, sortedTErrors, sortedAmbiguityErrors, sortedAmbiguityCorrectedErrors, sortedPCErrors, sortedPAErrors, sortedPMuErrors, sortedPTErrors, sortedPAmbiguityErrors, sortedPAmbiguityCorrectedErrors]
	

#Call the function for all noise levels

simulationFolder = 'Results/generic_noise0/'
[sortedMu0, sortedCErrors0, sortedAErrors0, sortedMuErrors0, sortedTErrors0, sortedAmbiguityScores0, sortedAmbiguityCorrectedErrors0,
 sortedPCErrors0, sortedPAErrors0, sortedPMuErrors0, sortedPTErrors0, sortedPAmbiguityErrors0, sortedPAmbiguityCorrectedErrors0] = sortData(simulationFolder)
simulationFolder = 'Results/generic_noise0.02/'
[sortedMu2, sortedCErrors2, sortedAErrors2, sortedMuErrors2, sortedTErrors2, sortedAmbiguityScores2, sortedAmbiguityCorrectedErrors2,
 sortedPCErrors2, sortedPAErrors2, sortedPMuErrors2, sortedPTErrors2, sortedPAmbiguityErrors2, sortedPAmbiguityCorrectedErrors2] = sortData(simulationFolder)
simulationFolder = 'Results/generic_noise0.04/'
[sortedMu4, sortedCErrors4, sortedAErrors4, sortedMuErrors4, sortedTErrors4, sortedAmbiguityScores4, sortedAmbiguityCorrectedErrors4,
 sortedPCErrors4, sortedPAErrors4, sortedPMuErrors4, sortedPTErrors4, sortedPAmbiguityErrors4, sortedPAmbiguityCorrectedErrors4] = sortData(simulationFolder)
simulationFolder = 'Results/generic_noise0.06/'
[sortedMu6, sortedCErrors6, sortedAErrors6, sortedMuErrors6, sortedTErrors6, sortedAmbiguityScores6, sortedAmbiguityCorrectedErrors6,
 sortedPCErrors6, sortedPAErrors6, sortedPMuErrors6, sortedPTErrors6, sortedPAmbiguityErrors6, sortedPAmbiguityCorrectedErrors6] = sortData(simulationFolder)
simulationFolder = 'Results/generic_noise0.08/'
[sortedMu8, sortedCErrors8, sortedAErrors8, sortedMuErrors8, sortedTErrors8, sortedAmbiguityScores8, sortedAmbiguityCorrectedErrors8,
 sortedPCErrors8, sortedPAErrors8, sortedPMuErrors8, sortedPTErrors8, sortedPAmbiguityErrors8, sortedPAmbiguityCorrectedErrors8] = sortData(simulationFolder)
simulationFolder = 'Results/generic_noise0.1/'
[sortedMu10, sortedCErrors10, sortedAErrors10, sortedMuErrors10, sortedTErrors10, sortedAmbiguityScores10, sortedAmbiguityCorrectedErrors10,
 sortedPCErrors10, sortedPAErrors10, sortedPMuErrors10, sortedPTErrors10, sortedPAmbiguityErrors10, sortedPAmbiguityCorrectedErrors10] = sortData(simulationFolder)


simulationFolder = 'Results/generic_random/'
[sortedMuR, sortedCErrorsR, sortedAErrorsR, sortedMuErrorsR, sortedTErrorsR, sortedAmbiguityScoresR, sortedAmbiguityCorrectedErrorsR,
 sortedPCErrorsR, sortedPAErrorsR, sortedPMuErrorsR, sortedPTErrorsR, sortedPAmbiguityErrorsR, sortedPAmbiguityCorrectedErrorsR] = sortData(simulationFolder)


#Have a function that bins by mu

def binValues(orderedMu, sortedErrors):
	binnedMu = []
	binnedError = []
	binnedStdAbove = []
	binnedStdBelow = []
	binnedAmbiguityErrors = []
	
	binSize = 10
	boxCount = len(orderedMu) / binSize
	offset = 0
	for boxInd in range(0, boxCount):
		
		errorValues = [sortedErrors[i] for i in range(offset,offset+binSize)]
		binnedMu.append(boxInd)
		meanError = sum(errorValues) / float(len(errorValues))
		binnedError.append(meanError)
	
		
		#Compute all values above and below the mean
		above = []
		below = []
		
		
		for value in range(0, len(errorValues)):
			if errorValues[value] > meanError:
				above.append(errorValues[value])
			else:
				below.append(errorValues[value])
			

		binnedStdAbove.append(np.std(errorValues))
		binnedStdBelow.append(np.std(errorValues))
		
		#Make sure that the sd is also reported per bin
		
		
		offset += binSize
	return [binnedMu, binnedError, binnedStdAbove, binnedStdBelow]


#Bin all values that we are interested in
[binnedMu, binnedCErrors0, binnedStdAboveC0, binnedStdBelowC0] = binValues(sortedMu0, sortedCErrors0)
[binnedMu, binnedCErrors2, binnedStdAboveC2, binnedStdBelowC2]  = binValues(sortedMu2, sortedCErrors2)
[binnedMu, binnedCErrors4, binnedStdAboveC4, binnedStdBelowC4]  = binValues(sortedMu4, sortedCErrors4)
[binnedMu, binnedCErrors6, binnedStdAboveC6, binnedStdBelowC6]  = binValues(sortedMu6, sortedCErrors6)
[binnedMu, binnedCErrors8, binnedStdAboveC8, binnedStdBelowC8]  = binValues(sortedMu8, sortedCErrors8)
[binnedMu, binnedCErrors10, binnedStdAboveC10, binnedStdBelowC10]  = binValues(sortedMu10, sortedCErrors10)

[binnedMu, binnedAErrors0, binnedStdAboveA0, binnedStdBelowA0] = binValues(sortedMu0, sortedAErrors0)
[binnedMu, binnedAErrors2, binnedStdAboveA2, binnedStdBelowA2]  = binValues(sortedMu2, sortedAErrors2)
[binnedMu, binnedAErrors4, binnedStdAboveA4, binnedStdBelowA4]  = binValues(sortedMu4, sortedAErrors4)
[binnedMu, binnedAErrors6, binnedStdAboveA6, binnedStdBelowA6]  = binValues(sortedMu6, sortedAErrors6)
[binnedMu, binnedAErrors8, binnedStdAboveA8, binnedStdBelowA8]  = binValues(sortedMu8, sortedAErrors8)
[binnedMu, binnedAErrors10, binnedStdAboveA10, binnedStdBelowA10]  = binValues(sortedMu10, sortedAErrors10)

[binnedMu, binnedMuErrors0, binnedStdAboveMu0, binnedStdBelowMu0] = binValues(sortedMu0, sortedMuErrors0)
[binnedMu, binnedMuErrors2, binnedStdAboveMu2, binnedStdBelowMu2]  = binValues(sortedMu2, sortedMuErrors2)
[binnedMu, binnedMuErrors4, binnedStdAboveMu4, binnedStdBelowMu4]  = binValues(sortedMu4, sortedMuErrors4)
[binnedMu, binnedMuErrors6, binnedStdAboveMu6, binnedStdBelowMu6]  = binValues(sortedMu6, sortedMuErrors6)
[binnedMu, binnedMuErrors8, binnedStdAboveMu8, binnedStdBelowMu8]  = binValues(sortedMu8, sortedMuErrors8)
[binnedMu, binnedMuErrors10, binnedStdAboveMu10, binnedStdBelowMu10]  = binValues(sortedMu10, sortedMuErrors10)

[binnedMu, binnedTErrors0, binnedStdAboveT0, binnedStdBelowT0] = binValues(sortedMu0, sortedTErrors0)
[binnedMu, binnedTErrors2, binnedStdAboveT2, binnedStdBelowT2]  = binValues(sortedMu2, sortedTErrors2)
[binnedMu, binnedTErrors4, binnedStdAboveT4, binnedStdBelowT4]  = binValues(sortedMu4, sortedTErrors4)
[binnedMu, binnedTErrors6, binnedStdAboveT6, binnedStdBelowT6]  = binValues(sortedMu6, sortedTErrors6)
[binnedMu, binnedTErrors8, binnedStdAboveT8, binnedStdBelowT8]  = binValues(sortedMu8, sortedTErrors8)
[binnedMu, binnedTErrors10, binnedStdAboveT10, binnedStdBelowT10]  = binValues(sortedMu10, sortedTErrors10)

#Do the same binning but then for the completely random case




#Group the values
binnedCErrors = [binnedCErrors0, binnedCErrors2, binnedCErrors4, binnedCErrors6, binnedCErrors8, binnedCErrors10]
binnedCStdAbove = [binnedStdAboveC0, binnedStdAboveC2, binnedStdAboveC4, binnedStdAboveC6, binnedStdAboveC8, binnedStdAboveC10]
binnedCStdBelow = [binnedStdBelowC0, binnedStdBelowC2, binnedStdBelowC4, binnedStdBelowC6, binnedStdBelowC8, binnedStdBelowC10]

binnedAErrors = [binnedAErrors0, binnedAErrors2, binnedAErrors4, binnedAErrors6, binnedAErrors8, binnedAErrors10]
binnedAStdAbove = [binnedStdAboveA0, binnedStdAboveA2, binnedStdAboveA4, binnedStdAboveA6, binnedStdAboveA8, binnedStdAboveA10]
binnedAStdBelow = [binnedStdBelowA0, binnedStdBelowA2, binnedStdBelowA4, binnedStdBelowA6, binnedStdBelowA8, binnedStdBelowA10]

binnedMuErrors = [binnedMuErrors0, binnedMuErrors2, binnedMuErrors4, binnedMuErrors6, binnedMuErrors8, binnedMuErrors10]
binnedMuStdAbove = [binnedStdAboveMu0, binnedStdAboveMu2, binnedStdAboveMu4, binnedStdAboveMu6, binnedStdAboveMu8, binnedStdAboveMu10]
binnedMuStdBelow = [binnedStdBelowMu0, binnedStdBelowMu2, binnedStdBelowMu4, binnedStdBelowMu6, binnedStdBelowMu8, binnedStdBelowMu10]

binnedTErrors = [binnedTErrors0, binnedTErrors2, binnedTErrors4, binnedTErrors6, binnedTErrors8, binnedTErrors10]
binnedTStdAbove = [binnedStdAboveT0, binnedStdAboveT2, binnedStdAboveT4, binnedStdAboveT6, binnedStdAboveT8, binnedStdAboveT10]
binnedTStdBelow = [binnedStdBelowT0, binnedStdBelowT2, binnedStdBelowT4, binnedStdBelowT6, binnedStdBelowT8, binnedStdBelowT10]

#Plotting function

def plotBinnedValues(binnedMu, binnedErrors, binnedStdAbove, binnedStdBelow, ymin, ymax, title):

	#For every array in the binned errors,
	#plot a line in the figure
	plt.figure()
	ax = axes()
	colors = ["#3498db", "#e74c3c", "#900090"]

	
	offset1 = 0.1
	offset2 = 0.2
	offset3 = 0.3
	offset4 = 0.4
	offset5 = 0.5
	offset6 = 0.6
	
	binnedMu1 = [i-offset1 for i in binnedMu]
	binnedMu2 = [i-offset2 for i in binnedMu]
	binnedMu3 = [i-offset3 for i in binnedMu]
	binnedMu4 = [i-offset4 for i in binnedMu]
	binnedMu5 = [i-offset5 for i in binnedMu]
	binnedMu6 = [i-offset6 for i in binnedMu]
	
	binnedMuOffsets = [binnedMu1, binnedMu2, binnedMu3, binnedMu4, binnedMu5, binnedMu6]
	
	legendLines = []
	labels = []
	#Plot the error for the simulations
	
	
		
	
	for binnedError in range(0, len(binnedErrors)):
		print binnedErrors[binnedError]
		print binnedStdBelow[binnedError]
		print binnedStdAbove[binnedError]
		
		newStdBelow = []
		#Make sure that the stds are corrected
		for std in range(0, len(binnedStdBelow[binnedError])):
			newStd = binnedStdBelow[binnedError][std]
			if binnedErrors[binnedError][std]-newStd < 0:
				newStd = abs(0-binnedErrors[binnedError][std])
			newStdBelow.append(newStd)
		
		p = plt.errorbar(binnedMuOffsets[binnedError], binnedErrors[binnedError], yerr=[newStdBelow, binnedStdAbove[binnedError]], label='Simulations', linewidth=2)
		legendLines.append(p[0])
		#labels.append('Simulations')
	
	ax.set_xticklabels(['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100'], rotation=90)
	ax.set_xticks(range(0,10))
	ax.set_xlabel('Mu (% tumor fraction)')
	ax.set_ylabel('Error')
		
	#ax.set_title(title)
	ax.set_xlim(-1,10)
	ax.set_ylim(ymin,ymax)
	labels = ['Sd = 0', 'Sd = 0.02', 'Sd = 0.04', 'Sd = 0.06', 'Sd = 0.08', 'Sd = 0.1'] #Remove random if we wish to export the normal figures
	plt.legend(legendLines, labels, numpoints=1, fontsize='12')
	plt.tight_layout()

	#plt.show()
	plt.savefig(title + '_S1.svg')

plotBinnedValues(binnedMu, binnedCErrors, binnedCStdAbove, binnedCStdBelow, 0, 1.4, 'C')
plotBinnedValues(binnedMu, binnedAErrors, binnedAStdAbove, binnedAStdBelow, 0, 1.4, 'A')
plotBinnedValues(binnedMu, binnedMuErrors, binnedMuStdAbove, binnedMuStdBelow, 0, 0.7, 'Mu')
plotBinnedValues(binnedMu, binnedTErrors, binnedTStdAbove, binnedTStdBelow, -0.5, 11, 'T')
	



