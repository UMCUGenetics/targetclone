"""
	The goal of this script is to make a normal.snp.txt file for theta matching the input tumor snp file.
	
	
	Give the script a vcf file with the original information (of the normal sample) and the file with the tumor.snp data.

"""

import sys
import re


vcfFile = sys.argv[2]
tumorParsedFile = sys.argv[1]
outFile = sys.argv[3]

#First pre-pare the VCF file and get all the required data

vcfEntries = dict() #keep a dictionary for quick searching
with open(vcfFile, 'r') as vcf:
	
	for line in vcf:
		if re.match("^#.*", line): #skip header
			continue
		
		#Obtain the read frequency of the variant and reference alleles
		splitLine = line.split("\t")
		sampleInfoField = splitLine[7]
		
		
		splitInfo = sampleInfoField.split(";")
		refCount = 0
		varCount = 0
		for field in splitInfo:
			splitField = field.split("=")
			if splitField[0] == 'RO': #reference allele	
				refCount = splitField[1]
			if splitField[0] == 'AO':
				varCount = splitField[1]
				
				splitVarCount = varCount.split(",")
				if len(splitVarCount) > 1:
					
					varCount = 0
					for i in splitVarCount:
						varCount += int(i)
					varCount = str(varCount)
				
		chromosome = splitLine[0]
		pos = int(splitLine[1]) - 1 #the positions in the theta input file are -1 for some reason
		
		catChr = chromosome + "_" + str(pos) #concatenate chromosome and position to search through the dictionary easily later
		
		vcfEntries[catChr] = [refCount, varCount]
		
#for every line in the tumor parsed file

positionList = dict() #keep a dictionary of all the positions so that we can map back easily which ones we are looking for. 

with open(outFile, 'wb') as out:
	#write header for the outfile
	out.write("#Chrm\tPos\tRef_Allele\tMut_Allele")
	out.write("\n")
	with open(tumorParsedFile, 'r') as tumorFile:
		lineCount = 0	
		for line in tumorFile:
	
			if lineCount < 1: #header
				lineCount += 1
				continue
		
			splitLine = line.split("\t")
			
			chromosome = splitLine[0]
			position = splitLine[1] #for some reason this is +1 from the vcf file
			#get the ref and mut data from the vcf file
			
			catChr = chromosome + "_" + position
			
			
			if catChr in vcfEntries:
				alleleCounts = vcfEntries[catChr]
				refCount = alleleCounts[0]
				varCount = alleleCounts[1]
			else:
				refCount = '0'
				varCount = '0'
			
			#Write the data to a new file
			
			out.write(chromosome + "\t" + position + "\t" + refCount + "\t" + varCount + "\n")
		
		
		
		

