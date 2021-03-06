files = dict(
	
	referenceFile = '../TargetClone/InternalData/reference_T3209.txt',
	outputDir = '../Results/noise0.02/',
	segmentationFile = '../TargetClone/InternalData/pq_segmentation.txt',
	simulationProbabilityFile = 'simulationProbabilities.txt',
	targetCloneInstance = '../TargetClone/InternalData/targetClone.pkl'
)

general = dict(
	numberOfSVs = 36,
	numberOfArmChanges = 20,
	malignantWcLosses = 10,
	malignantArmLosses = 10,
	malignantSVGains = 20,
	cycleSVGains = 2,
	cycleArmGains = 3,
	cycleArmLosses = 8,
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