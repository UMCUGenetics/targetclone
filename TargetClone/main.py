#Main script, start TargetClone from here
#How to use: main.py inputFile outputDirectory
#The algorithm will output the C matrix, A matrix, a tree structure containing the tree and the mu values

#Import the scripts from the settings folder
import sys

#Import all relevant libraries
import pickle
from shutil import copyfile
import numpy as np
from editDistance import EditDistance

from mu import Mu
from sample import Sample
from c import C
from combinations import CMuCombination
from laf import LAF
from alleles import Alleles
from distance import EventDistances
from sampleParser import SampleParser
from tree import Graph
from segmentations import Segmentation
from run import TargetClone
import settings
import os.path
import os
from listForest import generateOutput

import time
start_time = time.time()


#Check in the settings if kmin and kmax have been updated since the last time the method was run. If true, we need to re-initialize TC. Otherwise,
#we load the pre-initialized method.
def initializeSettingsHistory(targetClone):
	filepath = settings.files['settingsHistory']
	#load the file and update kmin and kmax
	newSettings = dict(kmin=settings.general['kmin'], kmax=settings.general['kmax'])
	with open(filepath, 'wb') as handle:
		pickle.dump(newSettings, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
	#Also re-create the TargetClone instance (we can use this one in the simulation data)
	with open('InternalData/targetClone.pkl', 'wb') as handle:
		pickle.dump(targetClone, handle, protocol=pickle.HIGHEST_PROTOCOL)

#only initialize TC if the settings have been updated since last time, otherwise load the TC instance
def initializeIfNecessary(): 
	filepath = settings.files['settingsHistory']
	
	if not os.path.isfile(filepath):
		initializeSettingsHistory()
		
	with open(filepath, 'rb') as handle:
		historySettings = pickle.load(handle)
	
		
	update = False	
	if historySettings['kmin'] != settings.general['kmin']:
		update = True
	if historySettings['kmax'] != settings.general['kmax']:
		update = True

	if update == False:
		print "Initializing from history..."
		pkl_file = open('InternalData/targetClone.pkl', 'rb')
		targetClone = pickle.load(pkl_file)
		pkl_file.close()
	else:
		print "Settings updated, re-initializing TargetClone..."
		targetClone = TargetClone()
		#make a new pickle instance
		initializeSettingsHistory(targetClone)
		
	return targetClone

#Check if the settings and the user input have valid values
def checkInputAndSettings(inFile, outputDir):
	
	#Check if we can read the infile
	if not os.path.isfile(inFile):
		print "ERROR: input file does not exist"
		exit(1)
	
	#Check if the output directory exists
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
		
	#Check all the settings for validity
	#kmin and kmax
	if int(settings.general['kmin']) < 0 or int(settings.general['kmax']) < 0:
		print "ERROR: kmin and kmax cannot be smaller than 0"
		exit(1)
	if int(settings.general['kmin']) >= 20 or int(settings.general['kmax']) >= 20:
		print "ERROR: kmin and kmax cannot be larger than 20"
		exit(1)
	if int(settings.general['kmin']) > int(settings.general['kmax']):
		print "ERROR: kmin cannot be larger than kmax"
		exit(1)
	#Check values of mu
	if int(settings.general['muMin']) < 0 or int(settings.general['muMax']) < 0:
		print "ERROR: muMin and muMax cannot be smaller than 0"
		exit(1)
	if int(settings.general['muMin']) > 100 or int(settings.general['muMax']) > 100:
		print "ERROR: muMin and muMax cannot be larger than 100"
		exit(1)
	if int(settings.general['muMin']) > int(settings.general['muMax']):
		print "ERROR: muMin cannot be larger than muMax"
		exit(1)
		
	if int(settings.general['maximumIterations']) < 0:
		print "ERROR: maximum iterations cannot be negative"
		exit(1)
	if float(settings.general['snvLowerBound']) < 0:
		print "ERROR: snvLowerBound cannot be negative"
		exit(1)
	if int(settings.general['precursorPloidy']) < 0:
		print "ERROR: precursorPloidy cannot be negative"
		exit(1)
	if int(settings.general['precursorAlleleACount']) < 0:
		print "ERROR: precursorAlleleACount cannot be negative"
		exit(1)
	if int(settings.general['precursorAlleleBCount']) < 0:
		print "ERROR: precursorAlleleBCount cannot be negative"
		exit(1)
	
	if int(settings.general['precursorAlleleACount']) + int(settings.general['precursorAlleleBCount']) is not int(settings.general['precursorPloidy']):
		print "ERROR: ploidy of the precursor does not match allele balance"
		exit(1)
	
	#Penalty values
	
	if float(settings.fst['pCWeight']) < 0:
		print "ERROR: pCWeight cannot be negative"
		exit(1)
	if int(settings.fst['lohPosLowerBound']) < 0:
		print "ERROR: lohPosLowerBound cannot be negative"
		exit(1)
	if int(settings.fst['notLohPosLowerBound']) < 0:
		print "ERROR: notLohPosLowerBound cannot be negative"
		exit(1)
	if float(settings.fst['lafLohUpperBound']) < 0:
		print "ERROR: lafLohUpperBound cannot be negative"
		exit(1)
	if int(settings.fst['alternatingAFPosLowerBound']) < 0:
		print "ERROR: alternatingAFPosLowerBound cannot be negative"
		exit(1)
	if float(settings.fst['alternatingAFPenalty']) < 0:
		print "ERROR: alternatingAFPenalty cannot be negative"
		exit(1)
	if float(settings.fst['lossOfConfidentLohPenalty']) < 0:
		print "ERROR: lossOfConfidentLohPenalty cannot be negative"
		exit(1)
	if float(settings.fst['lossOfMediumConfidentLohPenalty']) < 0:
		print "ERROR: lossOfMediumConfidentLohPenalty cannot be negative"
		exit(1)
	if float(settings.fst['lossOfSNVPenalty']) < 0:
		print "ERROR: lossOfSNVPenalty cannot be negative"
		exit(1)
	
	#also check files in settings, do these exist?
	if not os.path.isfile(settings.files['segmentationFile']):
		print "ERROR: segmentation file does not exist"
		exit(1)
	if not os.path.isfile(settings.files['excludedSNVs']):
		print "ERROR: excludedSNVs file does not exist"
		exit(1)
	if not os.path.isfile(settings.files['settingsHistory']):
		print "ERROR: settingsHistory file does not exist, please re-generate this file manually in main.py"
		exit(1)

#Initialize TargetClone instance from file to save time	
targetClone = initializeIfNecessary()

#Input file
inFile = sys.argv[1]
outputDir = sys.argv[2]

#perform a check of all the inputs, is it ok?
checkInputAndSettings(inFile, outputDir)


#Input the sample names as well. Add these as a header. 
tmpSamples = SampleParser().parseSampleTxt(inFile)

#We use the same p/q segmentation throughout the entire program. The segmentation is a property of the measurements
segmentation = Segmentation()
segmentation.setSegmentationFromFile(settings.files['segmentationFile'])

#Make a basic tree where every tumor sample has the same precursor
measurementLength = len(tmpSamples[0].measurements.measurements)
somVarNum = len(tmpSamples[0].somaticVariants)

#Obtain the ploidy of the precursor
precursorPloidy = int(settings.general['precursorPloidy'])

precursorAlleleACount = int(settings.general['precursorAlleleACount'])
precursorAlleleBCount = int(settings.general['precursorAlleleBCount'])

#Check if the ploidy is correct 
totalPloidy = precursorAlleleACount + precursorAlleleBCount

precursorTumorFrequency = 100 #Check if the ploidy is different from 2 (or allele balance). If true, then the precursor is not a healthy cell and we have 100% tumor
if precursorAlleleACount != 1 or precursorAlleleBCount != 1 or precursorPloidy != 2:
	precursorTumorFrequency = 0 #In this case we have 100% tumor in the precursor. Only if the precursor it is normal it is 0% tumor. 

#Initialize the 'healthy' sample, this can now also be a precursor	
healthySample = Sample(None, None)
healthySample.C = [C([2,precursorPloidy])]*measurementLength
healthySample.A = [Alleles(precursorAlleleACount,precursorAlleleBCount)]*measurementLength
healthySample.Mu = [Mu(precursorTumorFrequency)]
#obtain the chromosome, start and end information from the other samples
healthySample.measurements = LAF([0.5]*measurementLength, tmpSamples[0].measurements.chromosomes, tmpSamples[0].measurements.starts, tmpSamples[0].measurements.ends) 
healthySample.somaticVariants = [0]*somVarNum
healthySample.somaticVariantsInd = tmpSamples[0].somaticVariantsInd
healthySample.setParent(None)
healthySample.name = 'Precursor' #do not call it healthy, it may also be a 4N precursor. 

#Make a dummy bestCMu for the healthy sample
eventDistances = targetClone.eventDistances
bestCMuHealthy = CMuCombination(C([2,precursorPloidy]), Mu(precursorTumorFrequency), eventDistances)

healthySample.bestCMu = [bestCMuHealthy]*measurementLength
healthySample.originalCMu = healthySample.bestCMu

#Define all samples
samples = [healthySample] + tmpSamples

#get all somatic variants in one binary numpy array
allSomaticVariants = np.zeros((somVarNum, len(samples)))
for sample in range(0, len(samples)):
	allSomaticVariants[:,sample] = samples[sample].somaticVariants
	
	#also set the segmentation while we're at it
	samples[sample].measurements.segmentation = segmentation
	
targetClone.segmentation = segmentation
#The run part will also draw the trees automatically. 
[eCMatrix, eAMatrix, eSamples, trees, iterationMu, iterationMessages] = targetClone.run(samples)

#Output the edge data in a list for every iteration

#for every tree, output:
#-name list
#-parent
#-weights
#-annotations
nameList = dict()
parentList = dict()
weightList = dict()
annotationList = dict()
muList = iterationMu

for treeInd in trees.keys():
	nameList[treeInd] = []
	parentList[treeInd] = []
	weightList[treeInd] = []
	annotationList[treeInd] = []
	tree = trees[treeInd]
	
	#get the edges
	for edgeInd in range(0, len(tree.edgeList)):
		edge = tree.edgeList[edgeInd]
		#edge is a number, link it to the right sample and the name
		weight = edge[0]
		weightList[treeInd].append(weight)
		node = edge[2]
		
		nameList[treeInd].append(node)
		parent = edge[1]
		#we need to find what the index of the parent is to obtain the edge annotations
		nodeInd = 0
		parentInd = 0
		for sampleInd in range(0, len(samples)):
			if samples[sampleInd].name == node:
				nodeInd = sampleInd
			if samples[sampleInd].name == parent:
				parentInd = sampleInd
		if type(node) == int: #this is the precursor
			nodeInd = len(samples)
		if type(parent) == int:
			parentInd = len(samples)
		parentList[treeInd].append(parent)
		
		annotationInd = (parentInd, nodeInd)
		annotation = tree.edgeAnnotations[annotationInd]
		annotationList[treeInd].append(annotation)
		
#Write all data to output files: cMatrix, aMatrix, inferred mu (can be part of the trees as well in an annotation!) and the trees.

#This annotation file is useful if we wish to see the data going in to the visualization, but is not necessary
# annotationFile = outputDir + '/annotationFile.txt'
# f = open(annotationFile, 'w')
# print >> f, "names = ", nameList
# print >> f, "parents = ", parentList
# print >> f, "annotations = ", annotationList
# print >> f, "weights = ", weightList
# print >> f, "mu = ", muList
# print >> f, "messages = ", iterationMessages
# 
# f.close()


cStringMatrix = np.empty(eCMatrix.shape, dtype=object)
for cloneInd in range(0, eCMatrix.shape[1]):
	for measurementInd in range(0,eCMatrix.shape[0]):
		if type(eCMatrix[measurementInd][cloneInd]) == float or eCMatrix[measurementInd][cloneInd] is None or np.isnan(eCMatrix[measurementInd][cloneInd]): #catch the nan value, this always changes because there are multiple nan types...
			cStringMatrix[measurementInd][cloneInd] = 'NA'
		else:
			cStringMatrix[measurementInd][cloneInd] = int(eCMatrix[measurementInd][cloneInd])
		
np.savetxt(outputDir + '/cMatrix.txt', cStringMatrix, fmt='%s', delimiter='\t') #write as a string to make sure that NA values can be present in the matrix. 

eAStringMatrix = np.empty(eAMatrix.shape, dtype=object)

for cloneInd in range(0, eAMatrix.shape[1]):
	for alleleInd in range(0,eAMatrix.shape[0]):
		if type(eAMatrix[alleleInd][cloneInd]) == float or eAMatrix[alleleInd][cloneInd] is None: #the allele is only a float when the value is nan
			eAStringMatrix[alleleInd][cloneInd] = 'NA'
		else:
			eAStringMatrix[alleleInd][cloneInd] = eAMatrix[alleleInd][cloneInd].getAllelesAsString()
		
np.savetxt(outputDir + 'aMatrix.txt', eAStringMatrix, fmt='%s', delimiter='\t')

#Also save the final mu and final graph to a file (such that the simulation errors can easily be computed)
f = open(outputDir + '/EstimatedTree.txt', 'w')
bestTree = trees[len(trees)-1]
f.write(str(bestTree.getGraph()))  # python will convert \n to os.linesep
f.close()

bestMu = iterationMu[len(iterationMu)-1]
np.savetxt(outputDir + '/EstimatedMu.txt', bestMu, fmt='%f', delimiter='\t')


#Generate the output trees visualization
generateOutput(nameList, parentList, annotationList, weightList, muList, iterationMessages, outputDir)



#Copy the warning icon of the trees to the output folder, this is safer than linking to a directory. 
copyfile(settings.files['warningIconLocation'], outputDir + '/warning_icon.png')

print("--- Finished in %s seconds ---" % (time.time() - start_time))