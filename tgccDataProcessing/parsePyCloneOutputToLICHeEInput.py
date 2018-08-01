"""
	Provide the directory with the PyClone output for 1 patient, and this file will be converted to a CP-based input for LICheE. 

"""

import sys

pycloneInFile = sys.argv[1]
licheeOutFile = sys.argv[2]

with open(licheeOutFile, 'w') as outF:

	#First read all sample names from the file, also get the first mutation ID
	
	previousId = ""
	initialChrom = ""
	initialPos = ""
	seenSamples = []
	licheeHeader = "#chr\tposition\tdescription\tNormal"
	with open(pycloneInFile, 'r') as inF:
		lineCount = 0
		previousId = ""
		for line in inF:
			line = line.strip()
			
			if lineCount < 1:
				lineCount = 1
				continue
			
			splitLine = line.split("\t")
			
			mutationId = splitLine[0]
			
			
			sampleId = splitLine[1]
			
			if sampleId in seenSamples:
				break
			
			seenSamples.append(sampleId)
			previousId = mutationId
			initialChrom = mutationId.split(":")[0]
			initialPos = mutationId.split(":")[1]
			licheeHeader = licheeHeader + "\t" + sampleId
	
	outF.write(licheeHeader)
	outF.write("\n")

	with open(pycloneInFile, 'r') as inF:
		lineCount = 0
		licheeLine = initialChrom + "\t" + initialPos + "\t" + previousId + "\t0.0\t"
		
		for line in inF:
			line = line.strip()
			
			if lineCount < 1:
				lineCount = 1
				continue
			
			splitLine = line.split("\t")
			
			mutationId = splitLine[0]
			sampleId = splitLine[1]
			cp = splitLine[3]
			
			chrom = mutationId.split(":")[0]
			pos = mutationId.split(":")[1]
			
			#when the mutation id changes, make a new line for lichee
			#as long as the mutation id is the same, write to the same line. 
			if mutationId == previousId:
				#The samples are always in the same order
				licheeLine += cp
				licheeLine += "\t" #maybe the trailing tab will give issues
			else:
				outF.write(licheeLine)
				outF.write("\n")
				#start a new line
				licheeLine = chrom + "\t" + pos + "\t" + mutationId + "\t0.0\t" + cp + "\t"
				
			previousId = mutationId	
			
			
			
			
	
	
