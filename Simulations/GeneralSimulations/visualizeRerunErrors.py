"""
	Show the average differences between the re-runs. We compute the error pairwise from each run to every other run of the same simulation data.
	Then we can take the average of all those differences. This can be done for C, A, mu and T. We can plot these averages in boxplots. 

	It may be better to do this in parallel. Then we have 100 jobs, each running on a different node. 

"""

import sys
from glob import glob
import computeTreeErrorOtherMetrics
from tree import Graph

mainDir = sys.argv[1] #where to read the reruns from
prefix = sys.argv[2]

def collectErrorsFromFile(file): #subdir
	text_file = open(file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines


cErrors = []
aErrors = []
muErrors = []
treeErrors = []

subdirs = glob(mainDir + "/" + prefix + "*")
for subdir in subdirs:
	
	if len(cErrors) > 100: #accidentally ran more than 100 iterations, we need just 100. 
		continue

	cErrorFiles = glob(subdir + "/cError.txt")
	aErrorFiles = glob(subdir + "/aError.txt")
	muErrorFiles = glob(subdir + "/muError.txt")
	
	#For the tree, read the real tree and the inferred tree to compute the ancestry swap error. 
	
	if len(cErrorFiles) < 1: #if there is no data, continue. Either all files or no files should be present
		continue
	
	cErrorFile = cErrorFiles[0]
	aErrorFile = aErrorFiles[0]
	muErrorFile = muErrorFiles[0]
	
	cError = collectErrorsFromFile(cErrorFile)[0]
	aError = collectErrorsFromFile(aErrorFile)[0]
	muError = collectErrorsFromFile(muErrorFile)[0]
	
	cErrors.append(cError)
	aErrors.append(aError)
	muErrors.append(muError)
	
	realTreeFile = glob(subdir + "/RealTrees_1.txt")[0]
	inferredTreeFile = glob(subdir + "/EstimatedTrees_1.txt")[0]
	
	stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile('RealTrees_1.txt', subdir)[0]
	tree = eval(stringDict)
	realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
	
	stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile("EstimatedTrees_1.txt", subdir)[0]
	tree = eval(stringDict)
	inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])

	#Compute the ancestry swap error
	[ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)

	summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
	ancestrySwapError = summedError / float(noOfSamplePairs)
	
	treeErrors.append(ancestrySwapError)


#Compute the pairwise error for each simulation
cErrorDifferences = []
aErrorDifferences = []
muErrorDifferences = []
treeErrorDifferences = []

cErrorDifferenceCounts = 0
aErrorDifferenceCounts = 0
muErrorDifferenceCounts = 0
treeErrorDifferenceCounts = 0

for simulationInd in range(0, len(cErrors)): #should have the same number of keys as the other dictionaries
	for simulationInd2 in range(simulationInd, len(cErrors)):
		#Within that simulation dataset, do a pairwise comparison	
		cDifference = abs(cErrors[simulationInd] - cErrors[simulationInd2])
		cErrorDifferences.append(cDifference)
		
		if cDifference > 0:
			cErrorDifferenceCounts += 1

		
		aDifference = abs(aErrors[simulationInd] - aErrors[simulationInd2])
		aErrorDifferences.append(aDifference)
		
		if aDifference > 0:
			aErrorDifferenceCounts += 1
		
		muDifference = abs(muErrors[simulationInd] - muErrors[simulationInd2])
		muErrorDifferences.append(muDifference)
		
		if muDifference > 0:
			muErrorDifferenceCounts += 1
	
		treeDifference = abs(treeErrors[simulationInd] - treeErrors[simulationInd2])
		treeErrorDifferences.append(treeDifference)
		
		if treeDifference > 0:
			treeErrorDifferenceCounts += 1
	

#Show the average difference

averageCDifference = sum(cErrorDifferences) / float(len(cErrorDifferences))
averageADifference = sum(aErrorDifferences) / float(len(aErrorDifferences))
averageMuDifference = sum(muErrorDifferences) / float(len(muErrorDifferences))
averageTreeDifference = sum(treeErrorDifferences) / float(len(treeErrorDifferences))

##ALTERNATIVE of just reporting how many repeats are different. 

#Not implemented

#Keep the average differences stored somewhere (pkl?)

#then we need a collector script that goes through all these pkl files, and then combines everything into a figure

differenceFile = mainDir + "/" + prefix + "_differences.txt"

with open(differenceFile, 'w') as outF:
	
	outF.write(str(averageCDifference))
	outF.write("\n")
	outF.write(str(averageADifference))
	outF.write("\n")
	outF.write(str(averageMuDifference))
	outF.write("\n")
	outF.write(str(averageTreeDifference))
	
	
#A final script can then collect all of these difference files. 




