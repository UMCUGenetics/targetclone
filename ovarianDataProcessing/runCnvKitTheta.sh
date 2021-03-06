#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=20:00:00
#$ -e cnvKit_e
#$ -o cnvKit_o

#Script to run the copy number calling pipeline for ovarian data.

#1. Obtain the directory where all the BAM files and VCF files are stored per sample

folder="$1"
normalDir="$2" #directory with the reference samples
refFastaFile="$3"
thetaPath="$4"

#1.2 Obtain the fie locations of the normal files
normalBamFile=`ls "$normalDir"*.bam`
normalVcfFile=`ls "$normalDir"*_filtered.vcf`

#2. For each folder, run the pipeline (can be in parallel if necessary, but it seems very fast so far)

for d in "$folder"/*/ ; do
	
	bamFile=`ls "$d"*.bam`
	vcfFile=`ls "$d"*_filtered.vcf`
	
	echo "bam file: $bamFile"
	
	if [ ! -f "$vcfFile" ]; then
		echo "File not found!"
		continue
	fi
	
	#Extract the tumor sample file name	
	filename=$(basename -- "$bamFile")
	tumorFileName="${filename%.*}"
	
	echo "Tumor file: $tumorFileName"
	
	if [ "$tumorFileName" = "ERS312128" ] #skip the normal sample
	then
		continue
	fi
	
	if [ "$tumorFileName" != "ERS312132" ] #skip already processed files
	then
		continue
	fi
	
	#2.1 Run CNVKit to do the segmentation
	
	#Skip CNVKit, has already been done for all samples. Only run TheTA2

	cnvkit.py batch "$bamFile" --normal "$normalBamFile" \
    --fasta "$refFastaFile" -m wgs \
    --output-reference "$d"/cnvKit/ref.cnn --output-dir "$d/cnvKit/"
	
#	#2.2 Generate the tumor.snp files for theta
	cnsFile=`ls "$d"/cnvKit/*.cns`
	
	echo "Vcf file: $vcfFile"
	
	echo "ref: $d/cnvKit/ref.cnn"
	cnvkit.py export theta "$cnsFile" -r "$d"/cnvKit/ref.cnn -v "$vcfFile" #appears that it is not possible to specify output folder
	
	##2.3 Generate the normal.snp.txt files for theta
	tumorSnpFile=`ls *.tumor.snp_formatted.txt`
	echo "$tumorSnpFile"
	python generateThetaNormalAndTumorSnpFile.py "$tumorSnpFile" "$normalVcfFile" "$vcfFile" "$d"/cnvKit/normal.snp_formatted.txt "$d"/cnvKit/tumor.snp_formatted_corrected.txt
	
	#Move the temporary files to the right directory
	mv *.tumor_snp_formatted.txt "$d"/cnvKit/
	mv *.normal_snp_formatted.txt "$d"/cnvKit/
	mv *.interval_count "$d"/cnvKit/
	
	####RUNNING THETA IS NOT POSSIBLE ON THE HPC, THE MULTIPROCESSING LIBRARY OF PYTHON GIVES A UNICODE ERROR THAT I DON'T KNOW HOW TO SOLVE. SO AFTER GENERATING THE THETA FILES, THETA NEEDS TO BE RUN ON MAC WHERE IT DOES WORK
	#AFTER THAT COPY THE RESULTS BACK
	
	#2.3 Run TheTA
	#intervalCountFile=`ls $d/cnvKit/*.interval_count`
	#echo "interval file: $intervalCountFile"
	
	
	#"$thetaPath"/bin/RunTHetA "$intervalCountFile" --TUMOR_FILE "$d"/cnvKit/tumor.snp_formatted_corrected.txt --NORMAL_FILE "$d"/cnvKit/normal.snp_formatted.txt --BAF --FORCE -n 2 -d "$d"/cnvKit/
	##
	###Theta also does not seem to want to output to a specific directory (you shitty tools >:( ) so I move them after they have been created.
	#mv "$tumorFileName"* "$d"/output
	
	
done
