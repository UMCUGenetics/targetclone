#Script to generate the CNVKit VCF file for the TGCC TS data

#Get the ref and alt counts of all SNVs of the TGCC data


#1. Load the TGCC data (manually for now)

#2. Filter eval > 3

dataset = w3209 #switch for which tumor type dataset to use
tumorType = "T3209"

passInd = which(dataset$Eval > 3)
qualityData = dataset[passInd,]

#3. Obtain the SNP rows

snpInd = which(qualityData$SomVar != "Pos")
snpData = qualityData[snpInd,]

#4. Make a VCF file

sample = "PAT-275_6107EC2_070"

outFile = "~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T6107/EC70_test.vcf"

sink(outFile)
cat("id\tgene\ta\td\tmu_r\tmu_v\n")
for(row in 1:dim(somVars)[1]){
  id = paste0('s', row-1)
  gene = paste0(somVars$Ch[row], "_", somVars$Start[row])
  
  #Create the ref counts, separate by comma for each sample
  #Coverage is always the second column, the var counts the third. There are 11 columns per sample.
  
  #Also make sure that the right samples are used and not the excluded ones 
  
  refCoverage = ""
  altCoverage = ""
  
  startInd = 12 #columns with sample info start here
  nextCovInd = startInd + 1
  nextVarCovInd = startInd + 2
  sampleColNum = 11
  allColNames = colnames(somVars)
  for(i in startInd:dim(somVars)[2]){
    
    currentSample = allColNames[nextCovInd]
    
    splitSampleName = strsplit(currentSample,split='|', fixed=TRUE)[[1]][1]
    
    if(splitSampleName %in% excludedSamples){
      nextCovInd = nextCovInd + sampleColNum
      nextVarCovInd = nextVarCovInd + sampleColNum
      
      next
    }
    
    coverage = somVars[row,nextCovInd]
    varCoverage = somVars[row,nextVarCovInd]
    refCoverage = coverage - varCoverage
    
    if(i == startInd){
      allRefCounts = refCoverage
      allVarCounts = varCoverage
    }else{
      allRefCounts = paste(allRefCounts, refCoverage, sep=",")
      allVarCounts = paste(allVarCounts, varCoverage, sep=",")  
    }
    
    nextCovInd = nextCovInd + sampleColNum
    nextVarCovInd = nextVarCovInd + sampleColNum
    
    if(nextCovInd > dim(somVars)[2]){
      break
    }
    
  }
  #Write this line to the table for phylowgs
  line = paste(id, gene, allRefCounts, allVarCounts, '0.999', '0.499', sep="\t")
  cat(line)
  cat("\n")
  
}

sink()



