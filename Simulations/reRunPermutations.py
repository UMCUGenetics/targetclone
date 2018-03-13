#Given a folder

#Loop through the directories in this folder
#First remove the previous permutation data that will interfere 
#In every directory, re-run 10 permutations using runPermutation.py using the identifier
#This will re-generate the permutation data.

import sys
import os
import re

simulationFolder = sys.argv[1]

for subdir, dirs, files in os.walk(simulationFolder):
	if subdir == simulationFolder: #we are not interested in the root folder
		continue
	
	#check the files, and if one of these matchesa permutation type, we remove it.
	for file in files:
		#if any of these has to do with permutations, remove it from disk
		if re.match('pCError', file) or re.match('pAError', file) or re.match('pMuError', file) or re.match('pCError', file) or re.match('pTreeError', file) or re.match('pAmbiguityError', file) or re.match('pAmbiguityCorrectedError', file):
			f = open(subdir + '/' + file, 'w')
			f.close()