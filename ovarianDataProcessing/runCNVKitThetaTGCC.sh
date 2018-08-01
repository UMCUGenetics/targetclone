#Script to run the copy number calling pipeline for TGCC data.

#1. Obtain the directory where all the BAM files and VCF files are stored per sample

folder="$1"
normalDir="$2" #directory with the reference samples
targetFile="$3" #place to find the targets
refFastaFile="$4"

#1.2 Obtain the fie locations of the normal files
normalBamFile=`ls "$normalDir"*.bam`
normalVcfFile=`ls "$normalDir"*.vcf`

#2. For each folder, run the pipeline (can be in parallel if necessary, but it seems very fast so far)

for d in "$folder"/*/ ; do
	bamFile=`ls "$d"*.bam`
	vcfFile=`ls "$d"*.snp.vcf`
	
	#Extract the tumor sample file name	
	filename=$(basename -- "$bamFile")
	tumorFileName="${filename%.*}"
	
	echo "Tumor file: $tumorFileName"
	if [ "$tumorFileName" = "T618_NS30" ]
	then
		continue
	fi
	
	if [ "$tumorFileName" = "T618_NS48" ]
	then
		continue
	fi
	
	#2.1 Run CNVKit to do the segmentation
	
	#Skip CNVKit, has already been done for all samples. Only run TheTA2
#	
#	cnvkit.py batch "$bamFile" --normal "$normalBamFile" \
#    --targets "$targetFile" --fasta "$refFastaFile" -m amplicon \
#    --output-reference ref.cnn --output-dir "$d/output/"
#	
#	#2.2 Generate the tumor.snp files for theta
	cnsFile=`ls "$d"/output/*.cns`
	
	echo "Vcf file: $vcfFile"
	
	cnvkit.py export theta "$cnsFile" -r ref.cnn -v "$vcfFile" #appears that it is not possible to specify output folder
	
	##2.3 Generate the normal.snp.txt files for theta
	tumorSnpFile=`ls *.tumor.snp_formatted.txt`
	echo "$tumorSnpFile"
	python generateThetaNormalFile.py "$tumorSnpFile" "$normalVcfFile" "$d"/normal.snp_formatted.txt
	
	#
	##2.3 Run TheTA
	intervalCountFile=`ls *.interval_count`
	echo "interval file: $intervalCountFile"
	
	
	THetA/bin/RunTHetA "$intervalCountFile" --TUMOR_FILE "$tumorSnpFile" --NORMAL_FILE "$d"/normal.snp_formatted.txt --BAF --NUM_PROCESSES 2 --FORCE -n 2
	#
	##Theta also does not seem to want to output to a specific directory (you shitty tools >:( ) so I move them after they have been created.
	mv "$tumorFileName"* "$d"/output
	
done
