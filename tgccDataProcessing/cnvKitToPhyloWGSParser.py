"""
	Parse the output of cnvKit to a format that looks like Battenberg.

	X	Chrom	Start	End	x	Pval	x	x	major_cn	minor_cn	cp

"""
import sys
from glob import glob
import numpy as np

#Input a main directory

#search for folders with the right name from within these main folders
mainDir = sys.argv[1]
cp = sys.argv[2]
hg19CoordinateFile = sys.argv[3]

#1. Read the hg19 coordinates

coordinates = []
with open(hg19CoordinateFile, 'r') as f:
	lineCount = 0
	for line in f:
		line = line.strip()
		if lineCount < 1:
			lineCount = 1
			continue
		splitLine  = line.split("\t")
		
		coordinates.append([splitLine[1], splitLine[2], splitLine[3]])
		
hg19coordinates = np.array(coordinates)


subdirs = glob(mainDir + "/*/")

for subdir in subdirs:
	cnFiles = glob(subdir + "/output/*.cns.call")
	if len(cnFiles) < 1:
		continue
	
	cnFile = cnFiles[0]
	phylowgsOut = subdir + "/output/phylowgs.txt"
	print phylowgsOut
	with open(phylowgsOut, 'w') as outF:
		outF.write("x\tChrom\tStart\tEnd\tx\tPval\tx\tx\tmajor_cn\tminor_cn\tcp\n")
		with open(cnFile, 'r') as inF:
	
			lineCount = 0
			for line in inF:
				line = line.strip()
				if lineCount < 1:
					lineCount += 1
					continue
			
				splitLine = line.split("\t")
				
				chromFull = splitLine[0]
				splitChr = chromFull.split("chr")
				chrom = splitChr[1]
				
				start = splitLine[1]
				end = splitLine[2]
				
				#We don't need to check everything per se, could be faster
				for coordinate in hg19coordinates:
					if chrom == coordinate[0]:
						if start > coordinate[1]:
							#add a line for the starting segment
							outLine = "x\t" + chrom + "\t" + coordinate[1] + "\t" + start + "\tx\t1\tx\tx\t" + str(1) + "\t" + str(1) + "\t" + cp + "\n"
							outF.write(outLine)
						if end < coordinate[2]:
							#add ending segment
							outLine = "x\t" + chrom + "\t" + end + "\t" + coordinate[2] + "\tx\t1\tx\tx\t" + str(1) + "\t" + str(1) + "\t" + cp + "\n"
							outF.write(outLine)

					
				
				cn1 = int(splitLine[7])
				cn2 = int(splitLine[8])
				
				if cn1 > cn2:
					majorCn = cn1
					minorCn = cn2
				else:
					minorCn = cn1
					majorCn = cn2
				
				
				outLine = "x\t" + chrom + "\t" + splitLine[1] + "\t" + splitLine[2] + "\tx\t1\tx\tx\t" + str(majorCn) + "\t" + str(minorCn) + "\t" + cp + "\n"	
				outF.write(outLine)
				
				#Append to the file the regions that were not reported by CNVKit. Use the hg19 coordinates for that.
				#If the start of the segment is after the start of the chr, or the end is before the end of the chr, then add these lines to the file for that chromosome with a default normal cn
				
				
				
				

		
		
		

