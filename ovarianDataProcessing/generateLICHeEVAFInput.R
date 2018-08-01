#Script to generate LICHeE VAF-based input for TGCC data




#2. Filter eval > 3

dataset = W0618 #switch for which tumor type dataset to use
tumorType = "T618"

passInd = which(dataset$Eval > 3)
qualityData = dataset[passInd,]

#3. Obtain the SNV rows

somVarInd = which(qualityData$SomVar == "Pos")
somVarData = qualityData[somVarInd,]

#4. Get all SomVars
tumorTypeInd = which(somVarData$TumOri == tumorType)
somVars = somVarData[tumorTypeInd, ]

#5. Get the sample names in the right format

#> foo <- vector(mode="list", length=3)
#> names(foo) <- c("tic", "tac", "toe")
#> foo[[1]] <- 12; foo[[2]] <- 22; foo[[3]] <- 33

excludedSamplesT6107 = c('PA2-298_T6107-TE-Ad23_039|Coverage', 'PAT-358_T6107-PBL_069|Coverage', 'PAT-276_T92-6107_046|Coverage', 'PAT-358_T6107-NS_073|Coverage', 
                         'PA2-298_T6107-MET-Ad27_042|Coverage', 'PA2-298_T6107-MET-Ad26_041|Coverage', 'PA2-298_T6107-MET-Ad27_042|Coverage',
                         'PAT-360_T6107-MET-Ad25_020|Coverage')
excludedSamplesT3209 = c('PA2-296_T3209-MET-Ad21_037|Coverage', 'PA2-298_T3209-MET-Ad22_038|Coverage', 'PAT-253_3209_031|Coverage', 'PAT-275_3209TE1_072|Coverage',
                         'PAT-276_N3209_037|Coverage', 'PAT-358_T3209-PBL_070|Coverage', 'PAT-358_T3209-NS_074|Coverage', 'PAT-276_T95-3209_047|Coverage')
excludedSamplesT1382 = c('PAT-276_N1382_039|Coverage', 'PAT-358_T1382-PBL_072|Coverage', 'PAT-358_T1382-NS_076|Coverage', 'PA2-298_T1382-IEC-Li25_043|Coverage',
                         'PAT-275_T97-1382_049|Coverage', 'PA2-295_T1382-EC-0031_027|Coverage', 'PA2-295_T1382-MET-Ad06_028|Coverage', 'PA2-295_T1382-MET-Ad07_029|Coverage',
                         'PA2-296_T1382-MET-0037_031|Coverage', 'PA2-296_T1382-MET-Ad18_036|Coverage', 'PAT-369_T1382-TE-0030_023|Coverage')
excludedSamplesT618 = c('PAT-275_618CIS_051|Coverage', 'PAT-275_618Norm_052|Coverage', 'PAT-275_618Intratub_053|Coverage', 'PAT-276_N618_038|Coverage',
                        'PAT-358_T0618-NAP_071|Coverage')

excludedSamples = excludedSamplesT618

startInd = 12 #columns with sample info start here
nextCovInd = startInd + 1
nextVarCovInd = startInd + 2
sampleColNum = 11
allColNames = colnames(somVars)
samples = c()
for(i in startInd:dim(somVars)[2]){
  
  currentSample = allColNames[nextCovInd]
  if(currentSample %in% excludedSamples){
    nextCovInd = nextCovInd + sampleColNum
    nextVarCovInd = nextVarCovInd + sampleColNum
    next
  }
  if(is.na(currentSample)){
    break
  }
  samples = c(samples, currentSample)
  nextCovInd = nextCovInd + sampleColNum
  nextVarCovInd = nextVarCovInd + sampleColNum
}

sampleMap <- vector(mode="list", length=length(samples))
names(sampleMap) <- samples
#Make the map for the specific tumors

#For T6107
sampleMap[[1]] = "CIS30"
sampleMap[[2]] = "FCIS31"
sampleMap[[3]] = "YST40"
sampleMap[[4]] = "EC70"
sampleMap[[5]] = "EC71"
sampleMap[[6]] = "PBL36"
sampleMap[[7]] = "EC80"
sampleMap[[8]] = "TE81"
sampleMap[[9]] = "EC21"

#T3209
sampleMap[[1]] = "EB26"
sampleMap[[2]] = "CIS73"
sampleMap[[3]] = "TE74"
sampleMap[[4]] = "EC75"
sampleMap[[5]] = "TE82"
sampleMap[[6]] = "TE84"
sampleMap[[7]] = "EC85"
sampleMap[[8]] = "TE86"
sampleMap[[9]] = "EB87"
sampleMap[[10]] = "YST88"
sampleMap[[11]] = "EB22"
sampleMap[[12]] = "EB23"
sampleMap[[13]] = "EB24"
sampleMap[[14]] = "EB25"
sampleMap[[15]] = "CIS32"
sampleMap[[16]] = "EB33"
sampleMap[[17]] = "EC19"
sampleMap[[18]] = "YST20"

#T1382
sampleMap[[1]] = "MET30"
sampleMap[[2]] = "MET32"
sampleMap[[3]] = "MET33"
sampleMap[[4]] = "MET34"
sampleMap[[5]] = "MET35"
sampleMap[[6]] = "MET17"
sampleMap[[7]] = "MET18"
sampleMap[[8]] = "TE21"
sampleMap[[9]] = "EC22"
sampleMap[[10]] = "MET24"

#T618
sampleMap[[1]] = "FCIS27"
sampleMap[[2]] = "BCIS28"
sampleMap[[3]] = "CIS29"
sampleMap[[4]] = "NS30"
sampleMap[[5]] = "NS48"
sampleMap[[6]] = "NS75"


#dynamically generate outFile name
header = paste0("#chr\tposition\tdescription\tNormal")
for(sample in sampleMap){
  header = paste0(header, "\t", sample)
}

header = paste0(header, "\n")
#cat(header)




#It is better to first go through the samples, and then the rows to generate the files in one go. 
startInd = 12 #columns with sample info start here
nextCovInd = startInd + 1
nextVarCovInd = startInd + 2
sampleColNum = 11
allColNames = colnames(somVars)
varsToSkip = c()

#First go across the rows, then the columns

outFile = paste0("~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_Data/", tumorType, "/", tumorType, ".lichee.vaf.tsv")
sink(outFile)
cat(header)

for(row in 1:dim(somVars)[1]){

  #chromosome
  chromFull = somVars$Chrom[row]
  splitChr = strsplit(chromFull,split='chr', fixed=TRUE)
  chrom = splitChr[[1]][2]
  #position
  pos = somVars$Start[row]
  
  description = paste0(chrom, ":", pos)
  
  vafs = "0"
  startInd = 12 #columns with sample info start here
  nextCovInd = startInd + 1
  nextVarCovInd = startInd + 2
  sampleColNum = 11
  allColNames = colnames(somVars)
  
  #currentSample = allColNames[nextCovInd]
  #newSampleName = sampleMap[[currentSample]]
  
  
  for(i in startInd:dim(somVars)[2]){
    
    currentSample = allColNames[nextCovInd]
    newSampleName = sampleMap[[currentSample]]
    
    if(currentSample %in% names(sampleMap)){
      #Compute the VAF
      coverage = somVars[row,nextCovInd]
      varCoverage = somVars[row,nextVarCovInd]
      vaf = varCoverage / coverage
      if(is.na(vaf)){
        vaf = NaN
      }
      vafs = paste0(vafs, "\t", vaf)
      
      
    }

    nextCovInd = nextCovInd + sampleColNum
    nextVarCovInd = nextVarCovInd + sampleColNum
  }
  
  line = paste0(chrom, "\t", pos, "\t", description, "\t", vafs)
  #print(line)
  
  
  cat(line)
  cat("\n")
  
  
  
  
  
  
}

sink()


