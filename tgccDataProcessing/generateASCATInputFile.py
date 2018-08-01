"""
	Generate the LogR and BAF file for ASCAT.

"""

#Input: folder to search for cnr files from CNVKit for each sample
#Also files with the BAF for each sample (VAF)

#1. Read the positions from the BAF files
#2. Get the LogR

#4. Output all sample data in the same LogR and BAF files

import sys
import os
from glob import glob
import numpy as np

dataFolder = sys.argv[1]
logROutFilePrefix = sys.argv[2] #I will append tumor or germline to these automatically
bafOutFilePrefix = sys.argv[3]
germlineSample = sys.argv[4] #keep the germline sample separate so that it can be written to a different file. 

#First for the tumor samples
directories = glob(dataFolder + '/*/')

#Location of the bafFile
bafFile = glob(dataFolder + "/*_tumorBAF.txt")[0] ####THese should instead be the corrected files
germlineBafFile = glob(dataFolder + "/*_germlineBAF.txt")[0] 

#Collect all sample data
allLogR = []
allBAF = []
germlineLogR = []
chromosomes = []
positions = []
samples = [] #order in which the samples were read, for writing back to LogR file
skippedSNPs = dict()
positions = [] #the positions that remain after filtering
chromosomes = []
iteration = 0
for directory in directories:

	germline = False
	if directory == germlineSample: #if this is a germline sample, don't append it to the tumor LogR and tumor BAF
		germline = True
		
		
		
	#Get the cnr files from here
	
	subdir = directory + "/output"
	cnrFile = glob(subdir + "/*.cnr")[0]
	
	#Also get the allele frequencies
	
	
	

	currentLogRs = [] #make this numpy array later for quick overlap searches
	
	print cnrFile
	with open(cnrFile, "r") as f:
		
		lineCount = 0
		for line in f:
			if lineCount < 1:
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
			
			chromosome = splitLine[0]
			start = int(splitLine[1])
			end = int(splitLine[2])
			logR = splitLine[5]
			
			currentLogRs.append([chromosome, start, end, logR])
	
	
	currentLogRs = np.array(currentLogRs, dtype='object')
	
	
	
	#For each SNP in the tumor.snp_formatted, get these positions and the corresponding LogR from the cnr file
	sampleBAF = []
	sampleLogR = []
	with open(bafFile, 'r') as f:
		
		lineCount = 0
		header = []
		for line in f:
			
			if lineCount < 1:
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
	
			chromosome = splitLine[1]
			pos = splitLine[2]
			

			matchingChrInd = np.where(currentLogRs[:,0] == 'chr' + chromosome)
			subset = currentLogRs[matchingChrInd,:][0]
			
			overlappingSegmentsStart = subset[:,1] < int(pos)
			overlappingSegmentsEnd = subset[:,2] > int(pos)
			
			overlappingSegments = overlappingSegmentsStart * overlappingSegmentsEnd
			#There should be only one segment the SNP is in
			segment = subset[overlappingSegments, :]
			
			#Skip SNPs that do not have a LogR from CNVKit			
			#also remember which SNPs these are, because then we also need to remove them from the BAF file
			if len(segment) < 1:
				skippedSNPs[pos] = 0
				continue
			else: #store the LogR
				sampleLogR.append(segment[0][3])
				
				if iteration == 0: #also save the positional information
					chromosomes.append(chromosome)
					positions.append(pos)
				
				
	
	iteration += 1
	if germline == False:
		allLogR.append(sampleLogR)
		sampleName = directory.split("/")[4] #Don't add the sample name if this is a germline sample, it will not be in the tumor files. 
		samples.append(sampleName)
	else:
		germlineLogR = sampleLogR

#After collecting the data write everything to files

#1. Make the LogR file

allLogRNp = np.transpose(np.array(allLogR))

positions = np.array(positions)
#chromosomes = np.array(chromosomes)

colNum = positions.shape[0]
rowNum = len(samples) + 3 #snpInd, chr, pos and all samples

logR = np.empty([colNum, rowNum], dtype='object')

snpInds = range(0,colNum)
logR[:,0] = [str(snpInd) for snpInd in snpInds]
logR[:,1] = chromosomes
logR[:,2] = positions

logR[:,3:logR.shape[0]] = allLogRNp

headerSamples = "\t".join(samples)
header = "snpInd\tchrs\tpos\t" + headerSamples + "\n"

with open(logROutFilePrefix + "_tumorLogR.txt", 'wb') as outFile:
	outFile.write(header)
	for rowInd in range(0, logR.shape[0]):
		
		line = "\t".join(list(logR[rowInd,:]))
		outFile.write(line)
		outFile.write("\n")
		

#2. Filter the BAF data to only regions without the SNPs.
savedLines = []
header = ""
snpInd = 0
with open(bafOutFilePrefix + "_tumorBAF2.txt", 'w') as out:
	with open(bafFile, 'r') as f:
		lineCount = 0
		for line in f:
			
			if lineCount < 1:
				header = line
				out.write(header)
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
			
			pos = splitLine[2]
			
			if pos not in skippedSNPs:
				splitLine[0] = str(snpInd)
				joinedLine = "\t".join(splitLine)
				out.write(joinedLine)
				snpInd += 1

#3. Make germline LogR (which is sample PBL_70. I think it should be separate from the rest for the pipeline to work)
#Each sample will have the same LogR from this one reference sample. It is just meant to provide a reference for each sample, where the normal is sufficient.
#The positions need to be filtered correctly like was done above.


germlineLogRMatrix = np.empty([colNum, rowNum], dtype='object')

snpInds = range(0,colNum)
germlineLogRMatrix[:,0] = [str(snpInd) for snpInd in snpInds]
germlineLogRMatrix[:,1] = chromosomes
germlineLogRMatrix[:,2] = positions

for colInd in range(3, germlineLogRMatrix.shape[1]): #fill each column with the germline LogR
	
	germlineLogRMatrix[:,colInd] = germlineLogR
	

headerSamples = "\t".join(samples)
header = "snpInd\tchrs\tpos\t" + headerSamples + "\n"
with open(logROutFilePrefix + "_germlineLogR.txt", 'wb') as outFile:
	outFile.write(header)
	for rowInd in range(0, germlineLogRMatrix.shape[0]):
		
		line = "\t".join(list(germlineLogRMatrix[rowInd,:]))
		outFile.write(line)
		outFile.write("\n")

#re-create the germline BAF file with the right SNP positions as well 


savedLines = []
header = ""
snpInd = 0
with open(bafOutFilePrefix + "_germlineBAF2.txt", 'w') as out:
	with open(germlineBafFile, 'r') as f:
		lineCount = 0
		for line in f:
			
			if lineCount < 1:
				header = line
				out.write(header)
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
			
			pos = splitLine[2]
			
			if pos not in skippedSNPs:
				splitLine[0] = str(snpInd)
				joinedLine = "\t".join(splitLine)
				out.write(joinedLine)
				snpInd += 1
