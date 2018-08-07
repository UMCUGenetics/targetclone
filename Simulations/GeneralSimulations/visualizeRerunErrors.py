"""
	Show the average differences between the re-runs. We compute the error pairwise from each run to every other run of the same simulation data.
	Then we can take the average of all those differences. This can be done for C, A, mu and T. We can plot these averages in boxplots. 

	It may be better to do this in parallel. Then we have 100 jobs, each running on a different node. 

"""

import sys
from glob import glob

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


cErrors = dict()
aErrors = dict()
muErrors = dict()

subdirs = glob(mainDir + "/" + prefix + "*")
for subdir in subdirs:
	
	#Check if the prefix is different. The errors for every prefix can be in a dictionary with the correct key.
	print subdir
	
	splitSubdir = subdir.split("_")
	prefix = splitSubdir[0] #the first element always matches the name of the intial run

	cErrorFiles = glob(subdir + "/cError.txt")
	aErrorFiles = glob(subdir + "/aError.txt")
	muErrorFiles = glob(subdir + "/muError.txt")
	
	#For the tree, read the real tree and the inferred tree to compute the ancestry swap error. 
	
	if len(cErrorFiles) < 1: #if there is no data, continue. Either all files or no files should be present
		continue
	
	cErrorFile = cErrorFiles[0]
	aErrorFile = aErrorFiles[0]
	muErrorFile = muErrorFiles[0]
	
	cError = collectErrorsFromFile(cErrorFile)
	aError = collectErrorsFromFile(aErrorFile)
	muError = collectErrorsFromFile(muErrorFile)
	
	if prefix not in cErrors:
		cErrors[prefix] = []
	if prefix not in aErrors:
		aErrors[prefix] = []
	if prefix not in muErrors:
		muErrors[prefix] = []
		
	cErrors[prefix].append(cError)
	aErrors[prefix].append(aError)
	muErrors[prefix].append(muError)
	
	
	# 
	# 	
	# 	
	# 
	# 
	# if re.match('RealTrees', file): #read the file and obtain the error
	# 	stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
	# 	tree = eval(stringDict)
	# 	realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
	# 	treeSizes.append(len(realTree.edgeList))
	# 
	# if re.match('EstimatedTrees', file): #read the file and obtain the error
	# 	stringDict = computeTreeErrorOtherMetrics.collectErrorsFromFile(file, subdir)[0]
	# 	tree = eval(stringDict)
	# 	inferredTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])
	# 	
	# 
	# #Compute the ancestry swap error
	# [ancestrySwapErrorAbsentInInferred, ancestrySwapErrorPresentInInferred, noOfSamplePairs] = computeTreeErrorOtherMetrics.computeAncestrySwapError(realTree, inferredTree)
	# 
	# summedError = (ancestrySwapErrorAbsentInInferred + ancestrySwapErrorPresentInInferred)
	# ancestrySwapErrors.append(summedError / float(noOfSamplePairs))	

#Compute the pairwise error for each simulation
print cErrors

#Store the errors per data type separately
allCErrors = []
allAErrors = []
allMuErrors = []

for simulationInd in range(0, len(cErrors)): #should have the same number of keys as the other dictionaries
	1+1
	#Within that simulation dataset, do a pairwise comparison
	
	
	
	
		
	




