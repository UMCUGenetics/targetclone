#Script to read data from existing simulations, and use this to re-run TargetClone.
import sys
sys.path.insert(0, '../TargetClone/')
from simulationDataParser import SimulationDataHandler
from simulationErrors import SimulationErrorHandler
import simulationSettings
from run import TargetClone
from mu import Mu
import pickle
from segmentations import Segmentation

#Read:
#- AF measurements
#- LAF measurements
#- Sample names

#Then run TargetClone

#Instantiate a Simulator object
#Read the true C, A, mu and T
#Compute the error of the simulation


#From a shell script, call TargetClone on one folder individually per node.

simulationFolder = sys.argv[1]

#Read the simulationData.pkl in this folder

pkl_file = open(simulationFolder + '/simulationData.pkl', 'rb')
simulationData = pickle.load(pkl_file)
pkl_file.close()

#Read the real data
cMatrix = simulationData.cMatrix
aMatrix = simulationData.aMatrix
realTree = simulationData.realTree

#Obtain the mu from the files
savedMu = []
with open(simulationFolder + '/RealMu_0.txt', "r") as muFile:
	
	for line in muFile:
		mu = float(line)
		savedMu.append(Mu(mu))

samples = simulationData.samples

#Run TargetClone on samples

pkl_file = open(simulationSettings.files['targetCloneInstance'], 'rb')
targetClone = pickle.load(pkl_file)
pkl_file.close()

#Set the segmentation
segmentation = Segmentation()
segmentation.setSegmentationFromFile(simulationSettings.files['segmentationFile'])

targetClone.segmentation = segmentation
[eCMatrix, eAMatrix, eSamples, trees, iterationMu, iterationMessages] = targetClone.run(samples)

simulationErrorHandler = SimulationErrorHandler()
allCScores = []
allAScores = []
allMuScores = []
allTreeScores = []

eCMatrixFloat = eCMatrix.astype(float)
cMatrixFloat = cMatrix.astype(float)
#cScore = self.computeCRMSE(cMatrixFloat, eCMatrixFloat)
cScore = simulationErrorHandler.computeCError(cMatrixFloat, eCMatrixFloat)

allCScores.append(cScore / cMatrixFloat.size)

#aScore = self.computeARMSE(aMatrix, eAMatrix)
aData = simulationErrorHandler.computeAError(aMatrix, eAMatrix)
aScore = aData[0] / float(aMatrix.size)

allAScores.append(aScore)
eAStringMatrix = aData[2]
aStringMatrix = aData[1]


muData = simulationErrorHandler.computeMuError(eSamples, savedMu)
muScore = muData[0]
allMuScores.append(muScore)
allRealMuT = muData[1]
allEMuT = muData[2]

treeScore = simulationErrorHandler.computeTreeError(trees, realTree)
allTreeScores.append(treeScore)
bestTree = trees[len(trees)-1]


print "all C errors: ", allCScores

print "all A errors: ", allAScores

print "all mu errors: ", allMuScores

print "all Tree errors: ", allTreeScores

f = open(simulationFolder + '/cError.txt', 'w')
f.write(str(allCScores[0]))  # python will convert \n to os.linesep
f.close()
f = open(simulationFolder + '/aError.txt', 'w')
f.write(str(allAScores[0]))  # python will convert \n to os.linesep
f.close()
f = open(simulationFolder + '/muError.txt', 'w')
f.write(str(allMuScores[0]))  # python will convert \n to os.linesep
f.close()
f = open(simulationFolder + '/treeError.txt', 'w')
f.write(str(allTreeScores[0])) # python will convert \n to os.linesep
f.close()