#The goal of this script is to given an input directory run the call for calling major and minor copy numbers for all samples using CNVKit.


folder="$1" #tumor folder


for d in "$folder"/*/ ; do
	
	#1. First get the cns file
	cnsFile=`ls "$d"/output/*.cns`
	
	#2. Also read the right vcf file
	vcfFile=`ls "$d"*.snp.vcf`
	
	#3. Get the purity from the theta output
	thetaOutFile=`ls "$d"/output/*.n2.results`
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
	echo "$cnsFile"
	
	echo "$vcfFile"
	
	cnvkit.py call "$cnsFile" -v "$vcfFile" --purity "$tumorMu" -y -m clonal -o "$cnsFile".call
	
done



