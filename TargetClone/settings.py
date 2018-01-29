general = dict(
	kmin = 1, #This is the minimum copy number value in C
	kmax = 6, #This is the maximum copy number value in C
	muMin = 0, #This is the minimum value of mu (percentages)
	muMax = 100, #This is the maximum value of mu (percentages)
	maximumIterations = 10, #The maximum number of iterations the model is allowed to make
	snvLowerBound = 0.04, #If an SNV has an AF below this frequency, we will consider it noise and not take it into account (considering it not-present). 
	precursorPloidy = 4, #The ploidy of the precursor in tree T1
	precursorAlleleACount = 2, #The distribution of parental (A and B) alleles in this precursor in T1
	precursorAlleleBCount = 2,
	mmSigma = 0.02, #This sigma was estimated from the noise in the reference samples of TGCC, and is used in the Gaussian Mixture Models to represent sequencing noise.
	treeUpdatingIterations = 50 #The number of times we start from the MSA and then randomly remove the 'worst' branches for which the same number of SNVs violate the ISA. 
)

#The file paths should be relative to the directory where main.py is!
files = dict(
	segmentationFile = 'InternalData/pq_segmentation.txt', #File containing the locations of the p and q arms on the reference genome (default hg19)
	excludedSNVs = '../Examples/excludedSNVs.txt', #File with chromosomes and positions of SNVs that should be excluded from analysis (for example, driver mutations, hotspots)
	settingsHistory = 'InternalData/settings_history.pickle', #Internal file in which the settings of the previous run are stored for quick initialization
	warningIconLocation = 'InternalData/warning_icon.png' #Location of the warning icon that will be shown below our trees
)

fst = dict(
	pCWeight = 10000, #The weight for the P(C|T) term. A higher value enforces stricter horizontal dependency, copy numbers are less likely to fluctuate.
	lohPosLowerBound = 10, #The number of consecutive (excluding NA) positions that should have a LAF lower than lafLohUpperBound to be considered LOH. 
	notLohPosLowerBound = 10, #The number of consecutive (excluding NA) positions that should have a LAF higher than lafLohUpperBound to be considered not LOH. 
	lafLohUpperBound = 0.3, #The highest value of a LAF to be considered potential LOH. 
	alternatingAFPosLowerBound = 10, #The number of consecutive (excluding NA) positions that should have a LAF higher than lafLohUpperBound in sample 1 and lower than 1-lafLohUpperBound in sample 2 vv. to be considered alternating alleles. 
	alternatingAFPenalty = float("inf"), #The distance penalty assigned to samples with alternating alleles. 
	lossOfConfidentLohPenalty = float("inf"), #The distance penalty assigned to samples where loss of LOH is confidently identified (support in measurements and inferred A matrix)
	lossOfMediumConfidentLohPenalty = 50, #The distance penalty assigned to samples where loss of LOH is not confidently identified (support in only inferred A matrix)
	lossOfSNVPenalty = 50 #The distance penalty assigned to samples when loss of an SNV is not tolerated (i.e. no evidence of loss in surrounding alleles).
)

trees = dict(
	precursor = False #Should the method attempt to add one precursor ('unsampled subclone') to the tree? This feature is in beta, use with caution!
	
)

