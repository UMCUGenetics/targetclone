"""
	Read the strelka files of all ovarian samples and write the VAFs of the SNVs to LICHeE input. 

"""

import sys
from glob import glob
import numpy as np
import re

#Input a main directory

#search for folders with the right name from within these main folders
mainDir = sys.argv[1]

subdirs = glob(mainDir + "/*/")
numSamples = len(subdirs) - 2 #one folder will be the results folder


#For LICHeE, we need a column with the SNVs per sample. Also, we need to make sure that variants are set to 0 when these are not shared. But how about missing data?
#I think I will remove the SNVs that are not shared between all samples, because we just cannot say if these are missed or not present.
#If not in the vcf, we skip the variant, if 0, then we measure 0. 

#How to format this correctly?
#Maybe first keep all positions, and then filter these positions that are not shared between all samples, so take the intersection.

variants = []
variantIndices = dict() #keep the pos and chr and an index of where these are in the matrix. Then we should check if it is already there, and fill in a value.
#Afterwards, we make it a numpy array, and for values that are empty, we use NA and filter for rows with NA. 
currentVariantIndex = 0
sampleCount = 0
for subdir in subdirs:
	vcfFiles = glob(subdir + "/*.snvs.filtered.vcf")
	if len(vcfFiles) < 1:
		continue
	
	vcfFile = vcfFiles[0]
	print "vcf file: ", vcfFile
	
	with open(vcfFile, 'r') as inF:
		
		for line in inF:
			
			if re.match("^#", line):
				continue
			
			splitLine = line.split("\t")
			
			chrom = splitLine[0]
			pos = splitLine[1]
			
			concatPos = chrom + "_" + pos
			if concatPos in variantIndices:
				rowIndex = variantIndices[concatPos]
			else: #initialize row with NA for each value
				variants.append([chrom, pos, concatPos, 0.0]) #0 for the normal sample
				variants[currentVariantIndex] += ["NaN"]*numSamples

				variantIndices[concatPos] = currentVariantIndex
				rowIndex = currentVariantIndex
				currentVariantIndex += 1
			#print "variant index: ", currentVariantIndex	
			
			ref = splitLine[3]
			alt = splitLine[4]
			
			formatField = splitLine[8]
			splitFormatField = formatField.split(":")
			
			auIndex = splitFormatField.index("AU")
			tuIndex = splitFormatField.index("TU")
			cuIndex = splitFormatField.index("CU")
			guIndex = splitFormatField.index("GU")
			
			tumorInfo = splitLine[10]
			splitTumorInfo = tumorInfo.split(":")
		
			#uuugh
			if ref == 'A':
				refCounts = splitTumorInfo[auIndex]
			if ref == 'T':
				refCounts = splitTumorInfo[tuIndex]
			if ref == 'C':
				refCounts = splitTumorInfo[cuIndex]
			if ref == 'G':
				refCounts = splitTumorInfo[guIndex]	
			
			refCount = refCounts.split(",")[0]
			
			if alt == 'A':
				altCounts = splitTumorInfo[auIndex]
			if alt == 'T':
				altCounts = splitTumorInfo[tuIndex]
			if alt == 'C':
				altCounts = splitTumorInfo[cuIndex]
			if alt == 'G':
				altCounts = splitTumorInfo[guIndex]			
			
			altCount = altCounts.split(",")[0]
			
			vaf = int(altCount) / float(int(altCount) + int(refCount))

			variants[rowIndex][sampleCount+4] = vaf
			
		sampleCount += 1 #next sample
	
	
variants = np.array(variants)
print variants

#Filter the variants
matches = (variants != 'NA').all(axis=1)
filteredVariants = variants[matches,:]

print filteredVariants

np.savetxt('../../Ovarian_data/Results/ovarian_lichee_vaf.txt', filteredVariants, fmt='%s', delimiter='\t', header='#chr\tposition\tdescription\tNormal\tI1\tI2\tII1\tII2\tIV1\tIV2\tIV3\tV1') 
