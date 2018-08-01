#Plotting copy number profiles inferred by TargetClone

setwd('/Users/mnieboer/Documents/Projects/TargetClone/Thesis2015/Literature/scripts')
load("Plotdata1013.RData")

plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 7), xaxt="n", yaxt="n",  cex=2, ylab="Copy number", xlab="",  cex.lab=1,main= "T6107_EC70" ) #

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 7, 1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes

#Plot the normal chromosome copy number

normalC = 2
segments = c(hg19.csize$offset, sum(as.numeric(hg19.csize$length)))
for(segment in 1:length(segments)){
  lines(c(segments[segment],segments[segment+1]),c(normalC,normalC), lwd=2)  
}

#Now add the TargetClone copy numbers to the plot

path="/Users/mnieboer/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/cnvKitTheta/chromFa/"
EC70_CN = read.table(file = paste0(path,'EC70_cn.txt'), sep = '\t', header = TRUE)

for(row in 1:dim(EC70_CN)[1]){
  points(EC70_CN$plotx[row], EC70_CN$cn[row], pch=16, cex=1.5, col="blue")
}



#Estimates by TheTA2


plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 7), xaxt="n", yaxt="n",  cex=2, ylab="Copy number", xlab="",  cex.lab=1,main= "T6107_EC70", ) #
mtext("Normal:0%Tumor:100%", cex=2)

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 7, 1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes

#Plot the normal chromosome copy number

normalC = 2
segments = c(hg19.csize$offset, sum(as.numeric(hg19.csize$length)))
for(segment in 1:length(segments)){
  lines(c(segments[segment]+5000000,segments[segment+1])-5000000,c(normalC-0.05,normalC-0.05), lwd=4)  
}


ec70_c = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 2, 2, 0, 2, 2, 0, 0)	
for(segment in 1:length(segments)){
  lines(c(segments[segment]+5000000,segments[segment+1]-5000000),c(ec70_c[segment]+0.05,ec70_c[segment]+0.05), lwd=4, col="blue")  
}

#Make the plot with LAF and sSNVs for this tumor to show that the THeTA2 results are incorrect. 


path="/Users/mnieboer/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/cnvKitTheta/chromFa/"
EC70_LAF = read.table(file = paste0(path,'EC70_laf.txt'), sep = '\t', header = TRUE)

plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 0.8), xaxt="n", yaxt="n",  cex=2, ylab="AF/LAF", xlab="",  cex.lab=1,main= "T6107_EC70" ) #

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 0.8, 0.1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes
legend("topright", c("SNP", "somatic SNV"), pch=c(16,1), col=c("black", "red"), cex = 0.5, pt.cex = 2.5)
for(row in 1:dim(EC70_LAF)[1]){
  
  if(EC70_LAF$somVar[row] == "Pos"){
    points(EC70_LAF$plotx[row], EC70_LAF$LAF[row], cex=0.75, lwd=3, col="red")
  }else{
    points(EC70_LAF$plotx[row], EC70_LAF$LAF[row], pch=16, cex=0.75)
  }
  
  
}


#Repeat but then for TE74 of T3209


plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 7), xaxt="n", yaxt="n",  cex=2, ylab="Copy number", xlab="",  cex.lab=1,) #

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 7, 1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes

#Plot the normal chromosome copy number

normalC = 2
segments = c(hg19.csize$offset, sum(as.numeric(hg19.csize$length)))
for(segment in 1:length(segments)){
  lines(c(segments[segment],segments[segment+1]),c(normalC,normalC), lwd=2)  
}

#Now add the TargetClone copy numbers to the plot

path="/Users/mnieboer/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/cnvKitTheta/chromFa/"
TE74_CN = read.table(file = paste0(path,'TE74_cn.txt'), sep = '\t', header = TRUE)

for(row in 1:dim(TE74_CN)[1]){
  points(TE74_CN$plotx[row], TE74_CN$cn[row], pch=16, cex=1.5, col="blue")
}

#theta estimates

plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 7), xaxt="n", yaxt="n",  cex=2, ylab="Copy number", xlab="",  cex.lab=1) #
mtext("Normal:0%Tumor:100%", cex=2)

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 7, 1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes

#Plot the normal chromosome copy number

normalC = 2
segments = c(hg19.csize$offset, sum(as.numeric(hg19.csize$length)))
for(segment in 1:length(segments)){
  lines(c(segments[segment]+5000000,segments[segment+1])-5000000,c(normalC-0.05,normalC-0.05), lwd=4)  
}


ec70_c = c(2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, NA, 2, 2, 2, 2)	
for(segment in 1:length(segments)){
  lines(c(segments[segment]+5000000,segments[segment+1]-5000000),c(ec70_c[segment]+0.05,ec70_c[segment]+0.05), lwd=4, col="blue")  
}


#Make the LAF/SNV plot


path="/Users/mnieboer/Documents/Projects/TargetClone/Code/TargetClone_PlosRevision/cnvKitTheta/chromFa/"
EC70_LAF = read.table(file = paste0(path,'TE74_laf.txt'), sep = '\t', header = TRUE)

plot(NA, type= "p", xlim=c(hg19.csize[pstart,2],sum(as.numeric(hg19.csize$length)[1:pend])),
     ylim=c(0, 0.8), xaxt="n", yaxt="n",  cex=2, ylab="AF/LAF", xlab="",  cex.lab=1) #

axis(1, at=(hg19.csize$length/2)+hg19.csize$offset, labels=rownames(hg19.csize), cex.axis=2,  las=2) #x-axis , las=2
axis(2, at=seq(0, 0.8, 0.1), las=2, cex.axis=1.5) 

abline(v=c(hg19.csize$offset, sum(as.numeric(hg19.csize$length))), col="grey") #vertical lines to separate chromosomes
legend("topright", c("SNP", "somatic SNV"), pch=c(16,1), col=c("black", "red"), cex = 0.5, pt.cex=2.5)
for(row in 1:dim(EC70_LAF)[1]){
  
  if(EC70_LAF$somVar[row] == "Pos"){
    points(EC70_LAF$plotx[row], EC70_LAF$LAF[row], cex=0.75, lwd=3, col="red")
  }else{
    points(EC70_LAF$plotx[row], EC70_LAF$LAF[row], pch=16, cex=0.75)
  }
  
  
}


