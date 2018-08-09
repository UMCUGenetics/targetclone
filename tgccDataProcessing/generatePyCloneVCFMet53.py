"""
	Take EC21 as input and collect these positions from the snp file for MET53.
	Then register the allele counts for all these SNVs. 

"""

import sys
import re

refVcfIn = sys.argv[1]
vcfIn = sys.argv[2]
vcfOut = sys.argv[3]


vcfPositions = dict()
with open(refVcfIn, 'r') as inF:
	
	for line in inF:
		
		if re.match("^#", line):
			continue
		
		splitLine = line.split("\t")
		
		chrom = splitLine[0]
		
		pos = splitLine[1]
		
		concat = chrom + "_" + pos
		vcfPositions[concat] = [chrom, pos]
		

with open(vcfOut, 'w') as outF:
	
	#write the specific header
	outF.write("##fileformat=VCFv4.1\n")
	outF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">\n')
	outF.write('##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of A alleles used in tiers 1,2">\n')
	outF.write('##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of T alleles used in tiers 1,2">\n')
	outF.write('##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of C alleles used in tiers 1,2">\n')
	outF.write('##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of G alleles used in tiers 1,2">\n')
	outF.write("##contig=<ID=chr1,assembly=hg19,length=249250621>\n")
	outF.write("##contig=<ID=chr10,assembly=hg19,length=135534747>\n")
	outF.write("##contig=<ID=chr11,assembly=hg19,length=135006516>\n")
	outF.write("##contig=<ID=chr12,assembly=hg19,length=133851895>\n")
	outF.write("##contig=<ID=chr13,assembly=hg19,length=115169878>\n")
	outF.write("##contig=<ID=chr14,assembly=hg19,length=107349540>\n")
	outF.write("##contig=<ID=chr15,assembly=hg19,length=102531392>\n")
	outF.write("##contig=<ID=chr16,assembly=hg19,length=90354753>\n")
	outF.write("##contig=<ID=chr17,assembly=hg19,length=81195210>\n")
	outF.write("##contig=<ID=chr18,assembly=hg19,length=78077248>\n")
	outF.write("##contig=<ID=chr19,assembly=hg19,length=59128983>\n")
	outF.write("##contig=<ID=chr2,assembly=hg19,length=243199373>\n")
	outF.write("##contig=<ID=chr20,assembly=hg19,length=63025520>\n")
	outF.write("##contig=<ID=chr21,assembly=hg19,length=48129895>\n")
	outF.write("##contig=<ID=chr22,assembly=hg19,length=51304566>\n")
	outF.write("##contig=<ID=chr3,assembly=hg19,length=198022430>\n")
	outF.write("##contig=<ID=chr4,assembly=hg19,length=191154276>\n")
	outF.write("##contig=<ID=chr5,assembly=hg19,length=180915260>\n")
	outF.write("##contig=<ID=chr6,assembly=hg19,length=171115067>\n")
	outF.write("##contig=<ID=chr7,assembly=hg19,length=159138663>\n")
	outF.write("##contig=<ID=chr8,assembly=hg19,length=146364022>\n")
	outF.write("##contig=<ID=chr9,assembly=hg19,length=141213431>\n")
	outF.write("##contig=<ID=chrM,assembly=hg19,length=16569>\n")
	outF.write("##contig=<ID=chrX,assembly=hg19,length=155270560>\n")
	outF.write("##contig=<ID=chrY,assembly=hg19,length=59373566>\n")
	outF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMET_YST53\n")
	with open(vcfIn, 'r') as inF:
		
		for line in inF:
			
			if re.match("^#", line):
				
				continue
			
			splitLine = line.split("\t")
			
			chrom = splitLine[0]
			splitChrom = chrom.split("chr")
			newChrom = splitChrom[1]
			pos = splitLine[1]
			
			
			#Check if this position is in the vcf positions of the other file
			#If true, keep this record
			
			concat = newChrom + "_" + pos
			
			if concat in vcfPositions:
				
				#Write to a slightly different format to match the rest of the files
				#AO=0;RO=7;DP=7	GT:AO:RO	0/1:0:7
				infoField = splitLine[7]
				splitInfo = infoField.split(";")
				refCount = 0
				varCount = 0
				for field in splitInfo:
					splitField = field.split("=")
					if splitField[0] == 'DP4': #allele counts
						alleleCounts = splitField[1]
					
						#Split the allele counts again to get ref-fw, ref-rv, alt-fw, alt-rv
						splitAlleleCounts = alleleCounts.split(",")
						
						refCount = int(splitAlleleCounts[0]) + int(splitAlleleCounts[1])
						varCount = int(splitAlleleCounts[2]) + int(splitAlleleCounts[3])
				
				
				dp = refCount + varCount
				id = splitLine[2]
				ref = splitLine[3]
				var = splitLine[4]
				
				au = "0,0"
				tu = "0,0"
				cu = "0,0"
				gu = "0,0"
				if ref == "A":
					au = str(refCount) + "," + str(refCount)
				if ref == "T":
					tu = str(refCount) + "," + str(refCount)
				if ref == "C":
					cu = str(refCount) + "," + str(refCount)
				if ref == "G":
					gu = str(refCount) + "," + str(refCount)
				
				if var == "A":
					au = str(varCount) + "," + str(varCount)
				if var == "T":
					tu = str(varCount) + "," + str(varCount)
				if var == "C":
					cu = str(varCount) + "," + str(varCount)
				if var == "G":
					gu = str(varCount) + "," + str(varCount)
				
				qual = "."
				filter = "PASS"
				sampleInfo = "SOMATIC"
				format = "DP:AU:CU:GU:TU"
				sampleFormat = str(dp) + ":" + str(au) + ":" + str(cu) + ":" + str(gu) + ":" + str(tu)
				
				
				newLine = newChrom + "\t" + pos + "\t" + id + "\t" + ref + "\t" + var + "\t" + qual + "\t" + filter + "\t" + sampleInfo + "\t" + format + "\t" + sampleFormat + "\n"
				
				outF.write(newLine)
				



