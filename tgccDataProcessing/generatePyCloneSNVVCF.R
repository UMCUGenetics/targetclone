#Script to generate VCF input for PhyloWGS for SNVs
#We will mimick strelka

#Thus, the var counts need to be in AU/TU/CU/GU
#Also the filter field needs to be empty. 



#1. Load the TGCC data

#2. Filter eval > 3

dataset = w6107 #switch for which tumor type dataset to use
tumorType = "T6107"

passInd = which(dataset$Eval > 3)
qualityData = dataset[passInd,]

#3. Obtain the SNV rows

somVarInd = which(qualityData$SomVar == "Pos")
somVarData = qualityData[somVarInd,]

#4. Get all SomVars
tumorTypeInd = which(somVarData$TumOri == tumorType)
somVars = somVarData[tumorTypeInd, ]

#4. Print the ref and alt counts to a PhyloWGS file VCF format strelka styled


#It is better to first go through the samples, and then the rows to generate the files in one go. 
startInd = 12 #columns with sample info start here
nextCovInd = startInd + 1
nextVarCovInd = startInd + 2
sampleColNum = 11
allColNames = colnames(somVars)
varsToSkip = c()
for(i in startInd:dim(somVars)[2]){
  
  currentSample = allColNames[nextCovInd]
  
  splitSampleName = strsplit(currentSample,split='|', fixed=TRUE)[[1]][1]
  
  if(splitSampleName %in% excludedSamples){
    nextCovInd = nextCovInd + sampleColNum
    nextVarCovInd = nextVarCovInd + sampleColNum
    
    next
  }
  
  #dynamically generate outFile name
  outFile = paste0("~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_Data/", splitSampleName, ".snv.pyclone.vcf")
  sink(outFile)
  cat("##fileformat=VCFv4.1")
  cat("\n")
  cat('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">')
  cat("\n")
  cat('##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of A alleles used in tiers 1,2">')
  cat("\n")
  cat('##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of T alleles used in tiers 1,2">')
  cat("\n")
  cat('##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of C alleles used in tiers 1,2">')
  cat("\n")
  cat('##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of G alleles used in tiers 1,2">')
  cat("\n")
  
  cat("##contig=<ID=chr1,assembly=hg19,length=249250621>")
  cat("\n")
  cat("##contig=<ID=chr10,assembly=hg19,length=135534747>")
  cat("\n")
  cat("##contig=<ID=chr11,assembly=hg19,length=135006516>")
  cat("\n")
  cat("##contig=<ID=chr12,assembly=hg19,length=133851895>")
  cat("\n")
  cat("##contig=<ID=chr13,assembly=hg19,length=115169878>")
  cat("\n")
  cat("##contig=<ID=chr14,assembly=hg19,length=107349540>")
  cat("\n")
  cat("##contig=<ID=chr15,assembly=hg19,length=102531392>")
  cat("\n")
  cat("##contig=<ID=chr16,assembly=hg19,length=90354753>")
  cat("\n")
  cat("##contig=<ID=chr17,assembly=hg19,length=81195210>")
  cat("\n")
  cat("##contig=<ID=chr18,assembly=hg19,length=78077248>")
  cat("\n")
  cat("##contig=<ID=chr19,assembly=hg19,length=59128983>")
  cat("\n")
  cat("##contig=<ID=chr2,assembly=hg19,length=243199373>")
  cat("\n")
  cat("##contig=<ID=chr20,assembly=hg19,length=63025520>")
  cat("\n")
  cat("##contig=<ID=chr21,assembly=hg19,length=48129895>")
  cat("\n")
  cat("##contig=<ID=chr22,assembly=hg19,length=51304566>")
  cat("\n")
  cat("##contig=<ID=chr3,assembly=hg19,length=198022430>")
  cat("\n")
  cat("##contig=<ID=chr4,assembly=hg19,length=191154276>")
  cat("\n")
  cat("##contig=<ID=chr5,assembly=hg19,length=180915260>")
  cat("\n")
  cat("##contig=<ID=chr6,assembly=hg19,length=171115067>")
  cat("\n")
  cat("##contig=<ID=chr7,assembly=hg19,length=159138663>")
  cat("\n")
  cat("##contig=<ID=chr8,assembly=hg19,length=146364022>")
  cat("\n")
  cat("##contig=<ID=chr9,assembly=hg19,length=141213431>")
  cat("\n")
  cat("##contig=<ID=chrM,assembly=hg19,length=16569>")
  cat("\n")
  cat("##contig=<ID=chrX,assembly=hg19,length=155270560>")
  cat("\n")
  cat("##contig=<ID=chrY,assembly=hg19,length=59373566>")
  cat("\n")
  
  header = paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", splitSampleName, "\n")
  cat(header)
  
  
  for(row in 1:dim(somVars)[1]){
    
    #chromosome
    chromFull = somVars$Chrom[row]
    splitChr = strsplit(chromFull,split='chr', fixed=TRUE)
    chrom = splitChr[[1]][2]
    #position
    pos = somVars$Start[row]
    ID = row
    
    ref = somVars$Ref[row]
    alt = somVars$Var[row]
    qual = '.'
    filter = 'PASS'
    
    format = 'DP:AU:CU:GU:TU'
    
    coverage = somVars[row,nextCovInd]
    varCoverage = somVars[row,nextVarCovInd]
    refCoverage = coverage - varCoverage
    
    
    
    
    #The vcf caller does not accept NA, so I can only use ,
    if(is.null(varCoverage) || length(varCoverage) < 1 ){
      varCoverage = '.'
    }else{
      if(is.na(varCoverage)){
        varCoverage = '.'
      }
    }
    if(is.null(refCoverage) || length(refCoverage) < 1 ){
      refCoverage = '.'
    }else{
      if(is.na(refCoverage)){
        refCoverage = '.'
      }
    }
    
    # if(paste0(chrom, ":", pos) %in% varsToSkip){ #keep variants that are 0 in some samples out of other samples as well to ensure thatall have the same SNVs. 
    #   next
    # }
    # 
    # if(varCoverage == 0 || varCoverage == '.'){
    #   if(refCoverage == 0 || refCoverage == '.'){
    #     varsToSkip = c(varsToSkip, paste0(chrom, ":", pos))
    #     next  
    #   }
    # }
    
    if(varCoverage == '.' || refCoverage == '.'){
      next
    }
    
    #sample = paste0('0/1:', varCoverage, ':', refCoverage)
    
    au = 0
    cu = 0
    gu = 0
    tu = 0
    if(ref == "A"){
      au = refCoverage
    }else if(ref == "C"){
      cu = refCoverage
    }else if(ref == "T"){
      tu = refCoverage
    }else if(ref == "G"){
      gu = refCoverage
    }
    
    if(alt == "A"){
      au = varCoverage
    }else if(alt == "C"){
      cu = varCoverage
    }else if(alt == "T"){
      tu = varCoverage
    }else if(alt == "G"){
      gu = varCoverage
    }
    
    
    
    dp = paste0(coverage, ":")
    au = paste0(au, ",", au, ":")
    cu = paste0(cu, ",", cu, ":")
    gu = paste0(gu, ",", gu, ":")
    tu = paste0(tu, ",", tu)
    sample = paste0(dp, au, cu, gu, tu)
    info = "SOMATIC"
    line = paste(chrom, pos, ID, ref, alt, qual, filter, info, format, sample, sep='\t')
    #Write this line to the vcf file
    
    cat(line)
    cat("\n")
  }
  sink()
  
  
  nextCovInd = nextCovInd + sampleColNum
  nextVarCovInd = nextVarCovInd + sampleColNum
  
}




