#Run ASCAT for the TGCC samples



library(ASCAT)
#load the data
dir = '~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T3209/'
ascatData = ascat.loadData(paste0(dir, 'T3209_tumorLogR.txt'), paste0(dir, 'T3209_tumorBAF2.txt'), Germline_LogR_file = paste0(dir, 'T3209_germlineLogR.txt'),
                           Germline_BAF_file = paste0(dir, 'T3209_germlineBAF2.txt'))

segmentation = ascat.asmultipcf(ascatData, ascat.gg = NULL, penalty = 5, wsample = NULL,
                                selectAlg = "exact", refine = TRUE)

ascat.plotSegmentedData(segmentation)
copyNumbers = ascat.runAscat(segmentation, pdfPlot=T)


#Try again for T6107
dir = '~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T6107/'
ascatData = ascat.loadData(paste0(dir, 'T6107_tumorLogR.txt'), paste0(dir, 'T6107_tumorBAF2.txt'), Germline_LogR_file = paste0(dir, 'T6107_germlineLogR.txt'),
                           Germline_BAF_file = paste0(dir, 'T6107_germlineBAF2.txt'))

segmentation = ascat.asmultipcf(ascatData, ascat.gg = NULL, penalty = 10, wsample = NULL,
                                selectAlg = "exact", refine = TRUE)

ascat.plotSegmentedData(segmentation)
copyNumbers = ascat.runAscat(segmentation, pdfPlot=T)

#T1382

dir = '~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T1382/'
ascatData = ascat.loadData(paste0(dir, 'T1382_tumorLogR.txt'), paste0(dir, 'T1382_tumorBAF2.txt'), Germline_LogR_file = paste0(dir, 'T1382_germlineLogR.txt'),
                           Germline_BAF_file = paste0(dir, 'T1382_germlineBAF2.txt'))

segmentation = ascat.asmultipcf(ascatData, ascat.gg = NULL, penalty = 10, wsample = NULL,
                                selectAlg = "exact", refine = TRUE)

ascat.plotSegmentedData(segmentation)
copyNumbers = ascat.runAscat(segmentation, pdfPlot=T)

#T618

dir = '~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T618/'
ascatData = ascat.loadData(paste0(dir, 'T618_tumorLogR.txt'), paste0(dir, 'T618_tumorBAF2.txt'), Germline_LogR_file = paste0(dir, 'T618_germlineLogR.txt'),
                           Germline_BAF_file = paste0(dir, 'T618_germlineBAF2.txt'))

segmentation = ascat.asmultipcf(ascatData, ascat.gg = NULL, penalty = 10, wsample = NULL,
                                selectAlg = "exact", refine = TRUE)


#Prepare the data for MEDICC, split into major and minor copy numbers

p3WithCoordinates = read.table('~/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/TGCC_data/T3209/T3209_tumorBAF2.txt')
segment = 0
copySegment = 0
previousValues = c()

#obtain the previous BAF values of the last segment to ensure that these are not of the same value as the first segment. 
for(i in 1:length(segmentation$Tumor_BAF_segmented)){
  previousValues = c(previousValues, segmentation$Tumor_BAF_segmented[[i]][dim(segmentation$Tumor_BAF_segmented[[1]])[1]])
}

#Get the previous segment. 
previousCopyNumbers = c()
for(i in 1:dim(copyNumbers$nA)[2]){
  previousCopyNumbers = c(previousCopyNumbers, copyNumbers$nA[1,i])
}

segmentChromosomes = c()
segmentBAF = c()
segmentPositions = c()
minorCNSegments = c()
majorCNSegments = c()
for(i in 1:dim(segmentation$Tumor_BAF_segmented[[1]])[1]){
  
  logR = c()
  #get all the values of this position from the individual sample lists
  for(j in 1:length(segmentation$Tumor_BAF_segmented)){
    logR = c(logR, segmentation$Tumor_BAF_segmented[[j]][i])
  }
  
  copyNumbersPosMinor = c()
  for(j in 1:dim(copyNumbers$nA)[2]){
    copyNumbersPosMinor = c(copyNumbersPosMinor, copyNumbers$nA[i,j])
  }
  copyNumbersPosMajor = c()
  for(j in 1:dim(copyNumbers$nB)[2]){
    copyNumbersPosMajor = c(copyNumbersPosMajor, copyNumbers$nB[i,j])
  }
  
  cnMatches = copyNumbersPos == previousCopyNumbers
  if(length(which(matches == FALSE)) > 0){
    copySegment = copySegment + 1
  }
  
  #Obtain the copy numbers for this sample
  #Compare to the previous line. If it is different, then it is a new segment. 
  matches = logR == previousValues
  
  if(length(which(matches == FALSE)) > 0){
    segment = segment + 1
    segmentBAF = rbind(segmentBAF, logR)
    snpPos = rownames(segmentation$Tumor_BAF_segmented[[1]])[i]
    
    #snpPos = i #does the index match with the sorting of the segmented BAF?
    #Obtain the chromosome that this segment is on
    chrInfo = p3WithCoordinates[which(p3WithCoordinates$V1 == snpPos),2] #this is still a factor
    posInfo = p3WithCoordinates[which(p3WithCoordinates$V1 == snpPos),3]
    
    chromosome= as.character(levels(chrInfo))[chrInfo]
    position= as.character(levels(posInfo))[posInfo]
    segmentPositions = c(segmentPositions, position)
    segmentChromosomes = c(segmentChromosomes, chromosome)
    
    minorCNSegments = rbind(minorCNSegments, copyNumbersPosMinor)
    majorCNSegments = rbind(majorCNSegments, copyNumbersPosMajor)
    
  }
  previousValues = logR
}

#rownames(segmentBAF) = segmentChromosomes
#2. Add sample names as column names
colnames(segmentBAF) = segmentation$samples

#1. Add column with positions (start of segment is sufficient)
segmentBAF = cbind(segmentBAF, segmentPositions)
segmentBAF= cbind(segmentBAF, segmentChromosomes)

#3. Write the matrix to a file (TC compatible)
#we could further generate this file by hand, needs chr, pos, reference, and then all the samples. 
#Later also merge with somatic variants. 

write.table(segmentBAF, 'P1_segmentBAF.txt', sep = "\t", quote=FALSE, row.names = FALSE)

#Now also output the major and minor copy numbers

#We need this information per chromosome
#For each chromosome, a new file is needed
#For the major and minor we also need separate files
#In each of these files, a sample is a new fasta line
#Within the fasta line, we have the cn of all the segments


#Issue: MEDICC is not working, possibly because the major copy numbers can be lower than the minor copy numbers?
#First check which one is lowest, assign this to the minor list. If tie, do not change. 
realMajorCNSegments = majorCNSegments
realMinorCNSegments = minorCNSegments
indToRemove = c()
for(row in 1:dim(minorCNSegments)[1]){
  
  #if a row has NA, skip that row
  if(length(which(is.na(minorCNSegments[row,]) == TRUE)) > 0){
    
    indToRemove = c(indToRemove, row)
    next
  }
  
  #Also skip a row if there is a CN > 4, we cannot use this for MEDICC.
  if(length(which(minorCNSegments[row,] > 4)) > 0){
    
    indToRemove = c(indToRemove, row)
    next
  }
  if(length(which(majorCNSegments[row,] > 4)) > 0){
    
    indToRemove = c(indToRemove, row)
    next
  }
  
  for(col in 1:dim(minorCNSegments)[2]){
    
    if(minorCNSegments[row,col] > majorCNSegments[row,col]){
      realMajorCNSegments[row,col] = minorCNSegments[row,col]
      realMinorCNSegments[row,col] = majorCNSegments[row,col]
    }else{
      realMajorCNSegments[row,col] = majorCNSegments[row,col]
      realMinorCNSegments[row,col] = minorCNSegments[row,col]
    }
    
  }
}





#First for the minor copy numbers
previousChr = -1
for(chr in 1:length(segmentChromosomes)){
  
  if(length(which((chr == indToRemove) == TRUE)) > 0){ #check if the row ishould be ignored
    print(chr)
    next
  }
  
  if(segmentChromosomes[chr] != previousChr){
    
    #A different chromosome, so a new file
    chrFile <- file(paste0("chr", segmentChromosomes[chr], "_minor.fasta"), "w")
    
    cat(paste0(">diploid"), file=chrFile)
    cat("\n", file=chrFile)
    cat(paste0(rep(1, length(currentCNs[,1])),sep="",collapse = ""), file=chrFile)
    cat("\n", file=chrFile)
    
    #Each column becomes a new line in the fasta file
    for(col in 1:dim(currentCNs)[2]){
      cat(paste0(">", ascatData$samples[col]), file=chrFile)
      cat("\n", file=chrFile)
      cat(paste(currentCNs[,col],sep="",collapse = ""), file=chrFile)
      cat("\n", file=chrFile)
    }
    
    close(chrFile)
    currentCNs = c()
  }
  
  currentCNs = rbind(currentCNs, realMinorCNSegments[chr,])
  
  
  previousChr = segmentChromosomes[chr]
}

#realMajorCNSegments[which(realMajorCNSegments > 4)] = 4 #set 5 to 4 for now, this is the maximum
#Repeat for the major copy numbers
previousChr = -1
for(chr in 1:length(segmentChromosomes)){
  
  if(length(which((chr == indToRemove) == TRUE)) > 0){ #check if the row ishould be ignored
    print(chr)
    next
  }
  
  if(segmentChromosomes[chr] != previousChr){
    print(paste0("writing: ", segmentChromosomes[chr]))
    #A different chromosome, so a new file
    chrFile <- file(paste0("chr", segmentChromosomes[chr], "_major.fasta"), "w")
    
    cat(paste0(">diploid"), file=chrFile)
    cat("\n", file=chrFile)
    cat(paste0(rep(1, length(currentCNs[,1])),sep="",collapse = ""), file=chrFile)
    cat("\n", file=chrFile)
    
    #Each column becomes a new line in the fasta file
    for(col in 1:dim(currentCNs)[2]){
      cat(paste0(">", ascatData$samples[col]), file=chrFile)
      cat("\n", file=chrFile)
      cat(paste(currentCNs[,col],sep="",collapse = ""), file=chrFile)
      cat("\n", file=chrFile)
    }
    
    close(chrFile)
    currentCNs = c()
  }
  
  currentCNs = rbind(currentCNs, realMajorCNSegments[chr,])
  
  
  previousChr = segmentChromosomes[chr]
}  

