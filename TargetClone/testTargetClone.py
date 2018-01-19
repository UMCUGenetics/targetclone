import os
import sys
import re
import numpy as np
from tree import Graph
from alleles import Alleles
from simulationErrors import SimulationErrorHandler

## Script to test if the output of TargetClone matches the ground truth
# use as: python testTargetClone.py input.txt outputDir/ groundTruthDir/

#load main.py, this will run TargetClone
command = "python main.py " + sys.argv[1] + " " + sys.argv[2]
command = command.strip()


os.system(command)

#Read output files back into matrices so that we can compare

eCMatrix = np.loadtxt(sys.argv[2] + '/cMatrix.txt', dtype=int)
eAMatrix = np.loadtxt(sys.argv[2] + '/aMatrix.txt', dtype=str)
eMu = np.loadtxt(sys.argv[2] + '/EstimatedMu.txt', dtype=float)

#Obtain the tree
text_file = open(sys.argv[2] + '/EstimatedTree.txt', "r")
lines = text_file.read()
stringDict = []
	
for line in lines.split("\n"):
	if line != "":
		stringDict.append(line)

text_file.close()


tree = eval(stringDict[0])
eTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])

#Also read the real matrices with the ground truth

cMatrix = np.loadtxt(sys.argv[3] + '/RealC.txt', dtype=int)
aMatrix = np.loadtxt(sys.argv[3] + '/RealA.txt', dtype=str)
realMu = np.loadtxt(sys.argv[3] + '/RealMu.txt', dtype=float)

#Obtain the tree
text_file = open(sys.argv[3] + '/RealTree.txt', "r")
lines = text_file.read()
stringDict = []
	
for line in lines.split("\n"):
	if line != "":
		stringDict.append(line)

text_file.close()


tree = eval(stringDict[0])
realTree = Graph(tree['vertices'], set(tree['edges']), tree['edges'])


#Compute error rates (C, A, mu, T)

simulationErrorHandler = SimulationErrorHandler()
eCMatrixFloat = eCMatrix.astype(float)
cMatrixFloat = cMatrix.astype(float)
cScore = simulationErrorHandler.computeCError(cMatrixFloat, eCMatrixFloat)

print "C error: ", cScore / cMatrixFloat.size

#The A matrices need to be converted to A object matrices.
aObjMatrix = np.empty(aMatrix.shape, dtype=object)
eAObjMatrix = np.empty(aMatrix.shape, dtype=object)
for row in range(0, aMatrix.shape[0]):
	for col in range(0, aMatrix.shape[1]):
		#generate allele object
		allele = aMatrix[row][col]
		AOccurrences = [m.start() for m in re.finditer('A', allele)]
		ACount = len(AOccurrences)
		BOccurrences = [m.start() for m in re.finditer('B', allele)]
		BCount = len(BOccurrences)
		
		alleleObj = Alleles(ACount, BCount)
		aObjMatrix[row][col] = alleleObj
		
		allele = eAMatrix[row][col]
		AOccurrences = [m.start() for m in re.finditer('A', allele)]
		ACount = len(AOccurrences)
		BOccurrences = [m.start() for m in re.finditer('B', allele)]
		BCount = len(BOccurrences)
		
		alleleObj = Alleles(ACount, BCount)
		eAObjMatrix[row][col] = alleleObj
				
aData = simulationErrorHandler.computeAError(aObjMatrix, eAObjMatrix)
print "A error: ", aData[0] / float(aMatrix.size)

muError = simulationErrorHandler.computeMuErrorFromVectors(realMu, eMu)
print "Mu error: ", muError

treeError = simulationErrorHandler.computeTreeError([eTree], realTree)
print "Tree error: ", treeError
