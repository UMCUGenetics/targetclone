"""
	Goal: take strelka file as input and filter out variants under a specified threshold.
	Filters: LowEVS and LowDepth
	
	For indels, an additional filter is HighDepth

"""
import sys
import re

inFile = sys.argv[1]
outFile = sys.argv[2] # fitlered variants

with open(outFile, 'w') as out:
	with open(inFile, 'r') as f:
		
		for line in f:
			
			
			if re.match("^#", line): #skip headers
				out.write(line)
				continue
			
			splitLine = line.split("\t")
			
			#FILTER is the 7th field
			
			filterField = splitLine[6]
			
			if filterField == "LowEVS" or filterField == "LowDepth" or filterField == "LowEVS;LowDepth":
				continue #skip these variants
			
			else:
				out.write(line)	
							
				
			
			
			
			
			
			
			
		
	
	
	
	
	
	
	
