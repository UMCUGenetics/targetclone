#!/bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=20:00:00
#$ -e bamProcessing_e
#$ -o bamProcessing_o

bamFolder="$1"
sampleName="$2"

#First merge all BAM files in the provided folder
#samtools merge "$bamFolder"/merged.bam "$bamFolder"/*.bam

#Then filter unmapped reads
#samtools view -bF 4 "$bamFolder"/merged.bam > "$bamFolder"/filtered.bam

#Replace read groups (picard)
#samtools addreplacerg -r ID:"$sampleName" -r LB:L1 -r SM:"$sampleName" -o "$bamFolder"/replaced.bam "$bamFolder"/filtered.bam

java -jar /hpc/cog_bioinf/ridder/users/mnieboer/Tools/picard.jar AddOrReplaceReadGroups RGID="$sampleName" RGLB=None RGPL=SOLID RGPU=None RGSM="$sampleName" I="$bamFolder"/filtered.bam O="$bamFolder"/replaced.bam

#Sort BAM file
samtools sort -o "$bamFolder"/"$sampleName".bam "$bamFolder"/replaced.bam

#Index the BAM file
samtools index "$bamFolder"/"$sampleName".bam

#Now in principle if this was succesful, the old BAMs and intermediates can be removed


