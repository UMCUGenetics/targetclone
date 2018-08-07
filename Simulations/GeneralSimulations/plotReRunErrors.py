"""
	Script to collect the rerun errors. First run visualizeRerunErrors.py (ont the best name) to create the difference files, is split across the cluster using computeReRunErrors.sh.
	Then use this script after that to collect all the files with the differences and to make a boxplot of those. 

"""


import sys
from glob import glob

mainDir = sys.argv[1] #where to read the reruns from


cDifferences = []
aDifferences = []
muDifferences = []
treeDifferences = []

differenceFiles = glob(mainDir + "/*_differences.txt")

for differenceFile in differenceFiles:

	#Read the errors from these files
	
	with open(differenceFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip()
			if lineCount == 0:
				cDifferences.append(line)
			if lineCount == 1:
				aDifferences.append(line)
			if lineCount == 2:
				muDifferences.append(line)
			if lineCount == 3:
				treeDifferences.append(line)
				
			lineCount += 1
			
	

#Make a boxplot of the results
print cDifferences
print treeDifferences
	
		
			
			
		

