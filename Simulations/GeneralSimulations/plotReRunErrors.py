"""
	Script to collect the rerun errors. First run visualizeRerunErrors.py (ont the best name) to create the difference files, is split across the cluster using computeReRunErrors.sh.
	Then use this script after that to collect all the files with the differences and to make a boxplot of those. 

"""


import sys
from glob import glob
from pylab import plot, figure, hold, boxplot, show, axes, xlim, savefig

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
				cDifferences.append(float(line))
			if lineCount == 1:
				aDifferences.append(float(line))
			if lineCount == 2:
				muDifferences.append(float(line))
			if lineCount == 3:
				treeDifferences.append(float(line))
				
			lineCount += 1
			
	

#Make a boxplot of the results

fig = figure()
ax = axes()
hold(True)

#Plot the boxplots
bp = boxplot(cDifferences, positions = [1], widths = 0.6)
bp = boxplot(aDifferences, positions = [3], widths = 0.6)
bp = boxplot(muDifferences, positions = [5], widths = 0.6)
bp = boxplot(treeDifferences, positions = [7], widths = 0.6)

xlim(0,8)
ax.set_xticklabels(['C', 'A', 'Mu', 'Trees'])
ax.set_xticks([1, 3, 5, 7])
ax.set_ylabel('Average difference')

#show()
savefig('rerun.svg')
			
			
		

