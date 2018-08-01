#Convert fasta with chr headers to headers without the chr annotation.


import sys
import re

inFile = sys.argv[1]
outFile = sys.argv[2]

with open(outFile, 'w') as outF:
	with open(inFile, 'r') as f:
		
		
		for line in f:
			
			if re.match("^>", line):
				
				
				#replace chr with nothing
				line = line.replace("chr", "")
				
			outF.write(line)	