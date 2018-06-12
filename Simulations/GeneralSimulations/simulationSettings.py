files = dict(

	outputDir = 'Results/snps_10000_noSNVs/',
	segmentationFile = '../../TargetClone/InternalData/pq_segmentation.txt',
	simulationProbabilityFile = "../../TargetClone/InternalData/lossGainProbabilityFile.txt",
	targetCloneInstance = '../../TargetClone/InternalData/targetClone.pkl'
)

general = dict(
	numberOfSNVs = 50,
	numberOfSNPs = 10000,
	kmin = 1,
	kmax = 6,
	cellCycles = 4,
	noiseLevel = 0.02, #was 0.03
	numberOfMixedClones = 4, #This is the number of subclones that we mix in on top of the major clone! The name of this setting needs to be better. 
	maximumMinorCloneFrequency = 50, #This value should change between simulations to see the effects
	minimumMinorCloneFrequency = 41, #To test between 0 and 10, we set this value to 1. A value of 0 makes no sense for this test.
	minimumSNVFrequency = 0.1, #We filter out every SNV below this frequency, as we assume that it comes from the minor clone and is thus contaminating. (Find optimal value for this frequency). This only works for the scripts that use the realVAF!
	randomTree = False
)
#
runType = dict(
	randomMeasurements = False, #Do we want random LAF and random SNVs to be assigned to each sample? 
	horizontalShuffle = False #Do we shuffle the LAF measurements randomly within a sample? (To test the influence of the horizontal dependency)
	
	
)
