"""
	Script to collect the rerun errors. First run visualizeRerunErrors.py (ont the best name) to create the difference files, is split across the cluster using computeReRunErrors.sh.
	Then use this script after that to collect all the files with the differences and to make a boxplot of those. 

"""


import sys
from glob import glob

mainDir = sys.argv[1] #where to read the reruns from


cErrors = []
aErrors = []
muErrors = []
treeErrors = []

differenceFiles = glob(mainDir + "/*_differences.txt")

for differenceFile in differenceFiles:

	#Read the errors from these files
	print differenceFile

