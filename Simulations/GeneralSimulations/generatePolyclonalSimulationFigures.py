#1. Load modules
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
import scipy as sp
import scipy.stats


plt.rcParams.update({'font.size': 14})

#2. Collect TC errors from the output files (this should be in one shared utility file)

def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines

def readData(dataFolder, noiseLevels, addition):
	
	groupedCErrors = dict()
	groupedAErrors = dict()
	groupedMuErrors = dict()
	groupedTreeErrors = dict()
	groupedAverageAncestrySwapErrors = dict()
	
	for noiseLevel in noiseLevels:
		simulationFolder = dataFolder + "_" + str(noiseLevel)
		
		#Read all the errors into one list for this noise level
		cErrors = []
		aErrors = []
		muErrors = []
		treeErrors = []
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
			
			summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
			averagedAncestrySwapError.append(summedError / float(noOfSamplePairs))
			
		#Gather the data per noise level
		groupedCErrors[noiseLevel] = cErrors
		groupedAErrors[noiseLevel] = aErrors
		groupedMuErrors[noiseLevel] = muErrors
		groupedTreeErrors[noiseLevel] = treeErrors
		groupedAverageAncestrySwapErrors[noiseLevel] = averagedAncestrySwapError
		
	#Return the grouped data
	return [groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAverageAncestrySwapErrors]


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
		
		groupedAboveStd.append(upper)
		groupedBelowStd.append(lower)
		
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

#Average errors
def averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestryErrors):

	averagedCErrors = averageData(groupedCErrors, 'C')
	averagedAErrors = averageData(groupedAErrors, 'A')
	averagedMuErrors = averageData(groupedMuErrors, 'Mu')
	
	averagedTreeErrors = averageData(groupedTreeErrors, 'T')
	averagedAncestryErrors = averageData(groupedAncestryErrors, 'C')

	return [averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAncestryErrors]

#3. Make line plots showing increase in errors across the increase in the number of subclones
# C, A, mu and T can be in one figure

#- Compute the average per number of subclones
#- Compute the standard deviations
#- Plot the lines (all in one figure, double axes)

def plotData(noiseLevels, errors, aboveStd, belowStd, labels, colorInd, lim, title):
	
	colors = ['#1f77b4', '#ff7f0e', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	
	#Add offsets to the x axis
	#actually these are not noise levels but # of subclones, chcange this later
	aOffset = 0
	bOffset = 0.03
	cOffset = 0.05
	dOffset = 0.07
	eOffset = 0.09
	
	offsetNoiseLevels = []
	
	offsetNoiseLevels.append([i - aOffset for i in noiseLevels])
	offsetNoiseLevels.append([i - bOffset for i in noiseLevels])
	offsetNoiseLevels.append([i - cOffset for i in noiseLevels])
	offsetNoiseLevels.append([i - dOffset for i in noiseLevels])
	offsetNoiseLevels.append([i - eOffset for i in noiseLevels])	
	
	offsets = [aOffset, bOffset, cOffset, dOffset, eOffset]
	
	
	contaminationLevels = ['0-10', '10-20', '20-30', '30-40', '40-50']
	
	#colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
	plt.figure()
	ax = plt.gca()

	legendLines = []
	labelsLegend = []
	colorInd = 0
	for contaminationInd in range(0, len(errors)):
		
		correctedBelowStd = []
		for std in range(0, len(belowStd[contaminationInd])):
			
			correctedBelowStd.append(errors[contaminationInd][std] - belowStd[contaminationInd][std])
		correctedAboveStd = []
		for std in range(0, len(aboveStd[contaminationInd])):
			
			correctedAboveStd.append(aboveStd[contaminationInd][std] - errors[contaminationInd][std])
		
		
		#Plot the error for the simulations
		#p = ax.errorbar(noiseLevels, errors[contaminationInd], yerr=[correctedBelowStd, correctedAboveStd], label='$E_C$', color=colors[colInd], linewidth=2)
		p = ax.errorbar(offsetNoiseLevels[contaminationInd], errors[contaminationInd], yerr=[correctedBelowStd, correctedAboveStd], label=contaminationLevels[contaminationInd], linewidth=2, color=colors[colorInd])
		legendLines.append(p[0])
		labelsLegend.append(contaminationLevels[contaminationInd])
		colorInd += 1
	#ax.legend(legendLines, labels, loc=2, numpoints=1)
	
	ax.set_ylim(lim[0],lim[1])
	ax.set_xlabel(r'Number of additional subclones in sample')
	ax.set_ylabel('Error')
	ax.set_xlim(0.5,4.5)
	ax.xaxis.set_ticks(np.arange(1, 5, 1))
	ax.legend(legendLines, labelsLegend, loc=2, numpoints=1)

	plt.tight_layout()
	
	#plt.show()
	plt.savefig(title)


def plotFigureOne(allCErrors, allAErrors, allMuErrors, allTreeErrors, allAncestryErrors, allCStdAbove, allAStdAbove, allMuStdAbove, allTreeStdAbove, allCStdBelow, allAStdBelow, allMuStdBelow, allTreeStdBelow, allAncestryStdAbove, allAncestryStdBelow):
	
	plotData(noiseLevels, allCErrors, allCStdAbove, allCStdBelow, ['Copy numbers'], 0, [0,1], 'Output/Polyclonality/polyclonal_C_all.svg')
	plotData(noiseLevels, allAErrors, allAStdAbove, allAStdBelow, ['Alleles'], 1, [0,1], 'Output/Polyclonality/polyclonal_A_all.svg')
	plotData(noiseLevels, allMuErrors, allMuStdAbove, allMuStdBelow, ['Tumor fraction'], 3, [0,0.6], 'Output/Polyclonality/polyclonal_Mu_all.svg')
	plotData(noiseLevels, allTreeErrors, allTreeStdAbove, allTreeStdBelow, ['Trees'], 4, [-1,17], 'Output/Polyclonality/polyclonal_T_all.svg')
	plotData(noiseLevels, allAncestryErrors, allAncestryStdAbove, allAncestryStdBelow, ['Trees'], 4, [0,0.3], 'Output/Polyclonality/polyclonal_Ancestry_all.svg')

	

def generatePolyclonalityFigure(simulationFolder, noiseLevels, contaminationLevels):
	
	allCErrors = []
	allAErrors = []
	allMuErrors = []
	allTreeErrors = []
	allAncestryErrors = []
	
	allCStdAbove = []
	allAStdAbove = []
	allMuStdAbove = []
	allTreeStdAbove = []
	allAncestryStdAbove = []
	
	allCStdBelow = []
	allAStdBelow = []
	allMuStdBelow = []
	allTreeStdBelow = []
	allAncestryStdBelow = []
	
	
	for contaminationLevel in contaminationLevels:
		dataFolder = simulationFolder + str(contaminationLevel)
		#F1. Read the data from all the simulation folders
		[groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestryErrors] = readData(dataFolder, noiseLevels, '')
		
		#F2. Average the errors per noise level
		[averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors, averagedAncestryErrors] = \
		averageErrorsPerNoiseLevel(groupedCErrors, groupedAErrors, groupedMuErrors, groupedTreeErrors, groupedAncestryErrors)
		
		allCErrors.append(averagedCErrors)
		allAErrors.append(averagedAErrors)
		allMuErrors.append(averagedMuErrors)
		allTreeErrors.append(averagedTreeErrors)
		allAncestryErrors.append(averagedAncestryErrors)
		
		#Obtain the standard deviation above and below the mean for the averages
		[groupedAboveStdC, groupedBelowStdC] = obtainStandardDeviations(groupedCErrors, averagedCErrors)
		[groupedAboveStdA, groupedBelowStdA] = obtainStandardDeviations(groupedAErrors, averagedAErrors)
		[groupedAboveStdMu, groupedBelowStdMu] = obtainStandardDeviations(groupedMuErrors, averagedMuErrors)
		[groupedAboveStdT, groupedBelowStdT] = obtainStandardDeviations(groupedTreeErrors, averagedTreeErrors)
		[groupedAboveStdAncestry, groupedBelowStdAncestry] = obtainStandardDeviations(groupedAncestryErrors, averagedAncestryErrors)
		
		allCStdAbove.append(groupedAboveStdC)
		allAStdAbove.append(groupedAboveStdA)
		allMuStdAbove.append(groupedAboveStdMu)
		allTreeStdAbove.append(groupedAboveStdT)
		allAncestryStdAbove.append(groupedAboveStdAncestry)
		
		allCStdBelow.append(groupedBelowStdC)
		allAStdBelow.append(groupedBelowStdA)
		allMuStdBelow.append(groupedBelowStdMu)
		allTreeStdBelow.append(groupedBelowStdT)
		allAncestryStdBelow.append(groupedBelowStdAncestry)

	#F3. Plot the error per noise level in one figure
	#plotFigureOne(averagedCErrors, averagedAErrors, averagedMuErrors, averagedTreeErrors,
	#			  groupedAboveStdC, groupedBelowStdC, groupedAboveStdA, groupedBelowStdA, groupedAboveStdMu, groupedBelowStdMu, groupedAboveStdT, groupedBelowStdT)

	plotFigureOne(allCErrors, allAErrors, allMuErrors, allTreeErrors, allAncestryErrors, allCStdAbove, allAStdAbove, allMuStdAbove, allTreeStdAbove, allCStdBelow, allAStdBelow, allMuStdBelow, allTreeStdBelow, allAncestryStdAbove, allAncestryStdBelow)
	
	return 0

#Make Figure 4A and S10

contaminationLevels = ['1_10', '11_20', '21_30', '31_40', '41_50']

simulationFolder = 'Results/mixedSubclones_' #Folder to look in for the results
#Noise levels can in this case be interpreted as the number of subclones in a sample
noiseLevels = [1, 2, 3, 4]

generatePolyclonalityFigure(simulationFolder, noiseLevels, contaminationLevels)



