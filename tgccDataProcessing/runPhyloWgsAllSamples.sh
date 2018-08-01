#First create the input for all samples in the given folder


folder="$1" #tumor folder
phylowgsPath="$2"

cnvPaths=""
snvPaths=""
snvTypes=""
sampleCounter=1
for d in "$folder"/*/ ; do
	
	#Get the sample name for the input file such that PyClone uses this sample name in the data
	filename=$(basename -- "$d")
	sampleName="${filename%.*}"
	if [ "$sampleName" = "PBL36" ] || [ "$sampleName" = "NAP38" ] || [ "$sampleName" = "PBL39" ] || [ "$sampleName" = "PBL37" ]; then #would be neater if all of these are provided by theuser
		echo "Excluding reference sample $sampleName"
		continue
	fi
	
	#1. Get the snv vcf file
	vcfFile=`ls "$d"/*.snv.phylowgs.vcf`
	
	if [ ! -f "$vcfFile" ]; then
		echo "File not found!"
		continue
	fi
	
	#2. Get the cnv file
	cnvFile=`ls "$d"/output/phylowgs.txt`
	
	#get the purity of the sample
	thetaOutFile=`ls "$d"/output/*.n2.results`
	
	if [ ! -f "$thetaOutFile" ]; then
		echo "File not found!"
		continue
	fi
	
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
	

	#3. Run the PhyloWGS input parser
	python "$phylowgsPath"/parser/parse_cnvs.py -f battenberg -c "$tumorMu" "$cnvFile" --cnv-output "$d"/output/cnvs.txt
	
	cnvPaths="$cnvPaths --cnvs sample$sampleCounter=$d/output/cnvs.txt"
	snvPaths="$snvPaths sample$sampleCounter=$vcfFile"
	snvTypes="$snvTypes --vcf-type sample$sampleCounter=strelka"
	sampleCounter=$((sampleCounter+1))
	
done

echo "$cnvPaths"
echo "$snvPaths"
echo "$snvTypes"

#create multi-sample input
#Somehow this command does not work, but when I manually copy paste it into the terminal it does work?? 
echo "python "$phylowgsPath"/parser/create_phylowgs_inputs.py "$cnvPaths" "$snvTypes" "$snvPaths" --output-cnvs "$folder"/cnv_data.txt --output-variants "$folder/"ssm_data.txt --output-params "$folder"/params.json"

#Then run PhyloWGS and do the data copying as well

echo "python "$phylowgsPath"/multievolve.py --num-chains 4 --params "$folder"/params.json --ssms "$folder"/ssm_data.txt --cnvs "$folder"/cnv_data.txt"
