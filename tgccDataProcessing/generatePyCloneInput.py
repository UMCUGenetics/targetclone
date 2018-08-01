"""
	Give script a vcf file in strelka format, and parse this to the PyClone tsv input format. 

"""

import sys
import re

vcfFile = sys.argv[1]
cnFile = sys.argv[2]
outFile = sys.argv[3]

with open(outFile, 'w') as outF:
	outF.write("mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n")
	with open(vcfFile, 'r') as inF:
		
		for line in inF:
			line = line.strip()
			
			if re.match("^#", line):
				continue
			
			splitLine = line.split("\t")
			
			chrom = splitLine[0]
			pos = splitLine[1]
			
			mutationId = chrom + ":" + pos
			
			ref = splitLine[3]
			alt = splitLine[4]
			
			formatField = splitLine[8]
			sampleData = splitLine[9]
			
			splitFormatField = formatField.split(":")
			splitSampleData = sampleData.split(":")
			
			auInd = splitFormatField.index("AU")
			tuInd = splitFormatField.index("TU")
			cuInd = splitFormatField.index("CU")
			guInd = splitFormatField.index("GU")
			
			if ref == "A":
				refCounts = splitSampleData[auInd]
			if ref == "T":
				refCounts = splitSampleData[tuInd]
			if ref == "C":
				refCounts = splitSampleData[cuInd]
			if ref == "G":
				refCounts = splitSampleData[guInd]
				
			if alt == "A":
				varCounts = splitSampleData[auInd]
			if alt == "T":
				varCounts = splitSampleData[tuInd]
			if alt == "C":
				varCounts = splitSampleData[cuInd]
			if alt == "G":
				varCounts = splitSampleData[guInd]
				
			refCount = refCounts.split(",")[0]
			varCount = varCounts.split(",")[0]
				
			#Overlap with the right CN, get the copy number (make sure that there is always overlap)
			majorCn = None
			minorCn = None
			with open(cnFile, 'r') as cnF:
				
				lineCount = 0
				for line in cnF:
					
					line = line.strip()
					if lineCount < 1:
						lineCount = 1
						continue
					
					#get the cn and see if the SNV is in this region
					
					#If not, then it is in a CNN region
					splitLine = line.split("\t")
					
					chromFull = splitLine[0]
					splitChr = chromFull.split("chr")
					cnChrom = splitChr[1]
					
					start = splitLine[1]
					end = splitLine[2]
					
					if chrom == cnChrom and start < pos and end > pos:
						cn1 = int(splitLine[7])
						cn2 = int(splitLine[8])
						
						if cn1 > cn2:
							majorCn = cn1
							minorCn = cn2
						else:
							minorCn = cn1
							majorCn = cn2
				
				#Check if there was a match, otherwise it is in a CNN region
				if majorCn is None or minorCn is None:
					if chrom == "X": #Assuming male, should be a param later
						majorCn = "1"
						minorCn = "0"
					else:
						majorCn = "1"
						minorCn = "1"
				
				#Write the PyClone input file
				
				outF.write(mutationId + "\t" + refCount + "\t" + varCount + "\t2\t" + str(minorCn) + "\t" + str(majorCn) + "\n")
					
			