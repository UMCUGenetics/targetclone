#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=2:00:00
#$ -e snpFilter_e
#$ -o snpFilter_o

vcfIn="$1" #provide just the name of the file, not the extension so that we can do the gzip
vcfOut="$2"

#First unzip the vcf
gunzip "$vcfIn".gz

#Then run the filter
vcfutils.pl varFilter -d 30 -Q 20 -D 100 "$vcfIn" > "$vcfOut"

#Gzip the vcf
gzip "$vcfIn"




