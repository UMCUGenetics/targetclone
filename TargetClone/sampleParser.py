import csv
import re
from sample import Sample

class SampleParser:
	
	#Read a txt file with all the LAF information, parse these to individual sample objects
	def parseSampleTxt(self, file):
		#Step 1. Read the file by row and remove the rows where no change (0 or 1) is observed from the reference sample.
		#ignore the lines that are som var
		interestingLines = []
		noSomVarLines = []
		checkFile = []
		header = []
		with open(file, "r") as inFile:
			array = []
			lineNum = 0
			for line in inFile:
				
				line = line.replace(",", ".") #make sure that commas are always dots
				#check the data in this line.
				#Split the line by tab
				line = line.strip('\r\n')
				splitLine = re.split("\t", line)
				if lineNum > 0: #skip the header
					
					
					if splitLine[2] == 'Y': #if this line is a somatic variant line, we keep it
						interestingLines.append(splitLine)
						checkFile.append("\t".join(splitLine))
						continue
					
					sampleData = splitLine[4:len(splitLine)] #get the measurement values of the non-reference samples
					snps = []
					
					for measurement in sampleData:
						
						if measurement != "NA" and splitLine[3] != "NA": #skip the NA positions
							
							if float(splitLine[3]) > 0.1 and float(splitLine[3]) < 0.9: #if the reference does not have a very low presence of the reference or variant allele
								
								
								snps.append(measurement)
								
					#annotatedLine = splitLine
					
					if len(snps) > 0: #there is at least one potential SNP on this line, let's keep it
						interestingLines.append(splitLine)
						noSomVarLines.append(splitLine)
						#convert to laf
						splitLine2 = []
						for item in splitLine:
							if re.match("^\d+?\.\d+?$", item) is not None and float(item) > 0.5:
								splitLine2.append(str(1-float(item)))
							else:
								splitLine2.append(item)
						checkFile.append("\t".join(splitLine2))
				else: #we have the header, set the sample names
					header = splitLine
				
				lineNum += 1
		
		#Use the interestingLines list to parse the rest of the information
		samples = []
		
		colNum = 0
		chromosomes = []
		starts = []
		somaticVariants = []
		for column in zip(*[line for line in interestingLines]):
			if colNum == 0: #chromosome names
				chromosomes = column
				colNum += 1
				continue
			if colNum == 1: #start positions
				starts = column
				colNum += 1
				continue
			if colNum == 2: #indication of som var yes/no
				somaticVariants = column
				colNum += 1
				continue
			if colNum == 3: #indication of the reference sample, we do not do anything else with this one now
				colNum += 1
				continue
			sample = Sample(None, None)
			sample.setDataFromRawMeasurements(header[colNum], column, chromosomes, starts, somaticVariants)
			samples.append(sample)
			colNum += 1
		
		return(samples)