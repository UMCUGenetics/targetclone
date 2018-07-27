"""
	Script to parse Strelka VCF output to PyClone txt input
	
	

"""

import sys
import re
import numpy as np

inFile = sys.argv[1]
cnFile = sys.argv[2] #the file with the cn segments
outFile = sys.argv[3]

segments = []
with open(cnFile, 'r') as cnF:
	
	for line in cnF:
		line = line.strip()
		splitLine = line.split("\t")
		
		if splitLine[5] == "II": #quck and dirty focus on one sample for now
			
			segments.append(splitLine)
			
segments = np.array(segments)	

with open(outFile, 'w') as outF:
	
	outF.write("mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n")
	with open(inFile, 'r') as inF:
		
		for line in inF:
			
			if re.match("^#", line): #skip all headers
				continue
			
			#To parse strelka output: 
			#get the ref and alt bases
			#Then read the REFU and ALTU counts
			
			splitLine = line.split("\t")
			
			chrom = splitLine[0]
			pos = splitLine[1]
			
			mutId = chrom + ":" + pos
			
			ref = splitLine[3]
			alt = splitLine[4]
			
			tumorInfoMeta = splitLine[8]
			tumorInfo = splitLine[10]
			splitTumorInfoMeta = tumorInfoMeta.split(":")
			
			
			auIndex = splitTumorInfoMeta.index("AU")
			tuIndex = splitTumorInfoMeta.index("TU")
			cuIndex = splitTumorInfoMeta.index("CU")
			guIndex = splitTumorInfoMeta.index("GU")
			
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
			#Now also match the major/minor copy numbers to this. (use the ones from ascat)
			splitChr = chrom.split("chr")[1]
			minorCn = None
			majorCn = None
			for segment in segments:
				
				if splitChr == segment[0] and pos > segment[1] and pos < segment[2]:
					
					#match on this chromosome segment
					a = int(segment[3])
					b = int(segment[4])
					
					if a > b:
						majorCn = a
						minorCn = b
					else:
						majorCn = b
						minorCn = a
					
					break
			
			if minorCn is not None or majorCn is not None: #if there is no cn for a segment, skip that SNV
				outF.write(mutId + "\t" + refCount + "\t" + altCount + "\t" + "2" + "\t" + str(minorCn) + "\t" + str(majorCn) + "\n")
					
				
			
