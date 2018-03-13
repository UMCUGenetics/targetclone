#Read the simlation data folder
#Read which mu there are
#Report on which one is missing

import os
import re
import numpy as np

def collectErrorsFromFile(file, subdir): #subdir
	text_file = open(subdir + '/' + file, "r")
	lines = text_file.read()
	floatLines = []
	
	for line in lines.split("\n"):
		if line != "":
			floatLines.append(float(line))
	
	text_file.close()
	return floatLines

dataFolder = '../Results/'
noiseLevels = [0.02]

allPresentMu = []

for noiseLevel in noiseLevels:
	simulationFolder = dataFolder + '/noise' + str(noiseLevel)
	
	for subdir, dirs, files in os.walk(simulationFolder):
		if subdir == simulationFolder: #we are not interested in the root folder
			continue
		for file in files:
			if re.match('RealMu', file): #read the file and obtain the error
				mu = collectErrorsFromFile(file, subdir)
				allPresentMu.append(mu[1])
				

expectedMu = range(0,101)
expectedMu = [i/100 for i in expectedMu]

print len(expectedMu)
print len(allPresentMu)

print np.sort(allPresentMu)

for mu in expectedMu:
	
	if mu not in allPresentMu:
		print "mu: ", mu, " is missing"
	
			
