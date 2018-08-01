#Generate VCF files that contain the reference and variant calls

#We need a separate VCF file for each sample

#The names of the vcf files need to correspond to the BAM file names 



#1. Get all SNPs of the tumor type

dataset = W0618 #switch for which tumor type dataset to use
tumorType = "T618"

passInd = which(dataset$Eval > 3)
qualityData = dataset[passInd,]

#3. Obtain the SNP rows

snpInd = which(qualityData$TumOri == "SNP")
snpData = qualityData[snpInd,]

#4. Output the information to a VCF file per sample

#VCF format:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SampleName

#AO and RO are allele counts and reference counts


#T3209 excluded samples
excludedSamples = c('PAT-276_N3209_037', 'PAT-358_T3209-PBL_070', 'PAT-276_T95-3209_047', 'PAT-253_3209_031', 'PAT-358_T3209-NS_074')
#T6107 excluded samples
excludedSamples = c('PAT-276_N6107_036', 'PAT-358_T6107-PBL_069', 'PAT-276_T92-6107_046', 'PAT-358_T6107-NS_073', 'PA2-298_T6107-MET-Ad27_042')
#T1382 excluded samples
excludedSamples = c('PAT-276_N1382_039', 'PAT-358_T1382-PBL_072', 'PAT-275_T97-1382_049', 'PAT-358_T1382-NS_076', 'PA2-298_T6107-TE-Ad23_03', 'PAT-360_T6107-MET-Ad25_020')
#T618 excluded samples
excludedSamples = c('PAT-275_618Norm_052', 'PAT-275_618Intratub_05', 'PAT-275_618CIS_051')


#I need a map for the sample recoding
#I don't know how the fuck I do that in R and I have no internet connection

#It is better to first go through the samples, and then the rows to generate the files in one go. 
startInd = 12 #columns with sample info start here
nextCovInd = startInd + 1
nextVarCovInd = startInd + 2
sampleColNum = 11
allColNames = colnames(snpData)
for(i in startInd:dim(snpData)[2]){
  
  currentSample = allColNames[nextCovInd]
  
  splitSampleName = strsplit(currentSample,split='|', fixed=TRUE)[[1]][1]
  
  if(splitSampleName %in% excludedSamples){
    nextCovInd = nextCovInd + sampleColNum
    nextVarCovInd = nextVarCovInd + sampleColNum
    
    next
  }
  
  #dynamically generate outFile name
  outFile = paste0("~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_Data/", splitSampleName, ".snp.vcf")
  sink(outFile)
  cat("##fileformat=VCFv4.1")
  cat("\n")
  cat("##INFO=<ID=AO,Number=A,Type=Integer,Description='Alternate allele observations'>")
  cat("\n")
  cat("##INFO=<ID=RO,Number=1,Type=Integer,Description='Reference allele observations'>")
  cat("\n")
  cat("##INFO=<ID=DP,Number=1,Type=Integer,Description='Read Depth'>")
  cat("\n")
  cat("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
  cat("\n")
  cat("##FORMAT=<ID=AO,Number=A,Type=Integer,Description='Alternate allele observation count'>")
  cat("\n")
  cat("##FORMAT=<ID=RO,Number=1,Type=Integer,Description='Reference allele observation count'>")
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
  
  
  for(row in 1:dim(snpData)[1]){
  
  
  #chromosome
  chrom = snpData$Chrom[row]
  #position
  pos = snpData$Start[row]
  ID = row
  #ref and alt can be . for theta?
  ref = snpData$Ref[row]
  alt = snpData$Var[row]
  qual = '.'
  filter = 'PASS'
  
  format = 'GT:AO:RO'
  
  
  coverage = snpData[row,nextCovInd]
  varCoverage = snpData[row,nextVarCovInd]
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
  
  
  sample = paste0('0/1:', varCoverage, ':', refCoverage)
  ao = paste0("AO=", varCoverage, ";")
  ro = paste0("RO=", refCoverage, ";")
  dp = paste0("DP=", coverage)
  info = paste0(ao, ro, dp)
  

  
  line = paste(chrom, pos, ID, ref, alt, qual, filter, info, format, sample, sep='\t')
  #Write this line to the vcf file
  
  cat(line)
  cat("\n")
  }
  sink()

  
  nextCovInd = nextCovInd + sampleColNum
  nextVarCovInd = nextVarCovInd + sampleColNum
    
}
  
  


