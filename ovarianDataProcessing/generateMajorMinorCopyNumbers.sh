#The goal of this script is to given an input directory run the call for calling major and minor copy numbers for all samples using CNVKit.

#Ovarian data specific script

folder="$1" #tumor folder

for d in "$folder"/*/ ; do
	
	#1. First get the cns file
	cnsFile=`ls "$d"/cnvKit/*.cns`
	
	echo "cns file: $cnsFile"
	
	#2. Also read the right vcf file
	vcfFile=`ls "$d"*_filtered.vcf`
	echo "vcf file: $vcfFile"
	
	#3. Get the purity from the theta output
	thetaOutFile=`ls "$d"/cnvKit/*.n2.results`
	echo "$thetaOutFile"
	lineCount=0
	while IFS='' read -r line || [[ -n "$line" ]]; do
		if [ "$lineCount" -lt 1 ]
		then
			lineCount=1
			continue
		fi
		echo "Text read from file: $line"
		mu=`awk -F "\t" '{print $2}' <<< "$line"`
		tumorMu=`awk -F "," '{print $2}' <<< "$mu"`
		echo "$tumorMu"
		#Split the line and get the mu information 
		
	done < "$thetaOutFile"
	
	#4. Then call cnvkit
	
	cnvkit.py call "$cnsFile" -v "$vcfFile" --purity "$tumorMu" -x female -m clonal -o "$cnsFile".call
	
	break
	
done



