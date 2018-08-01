


folder="$1" #tumor folder

pyCloneInputFiles=""
fileNum=0
for d in "$folder"/*/ ; do
	
	#1. Get the snv vcf file
	vcfFile=`ls "$d"/*.pyclone.vcf`
	
	if [ ! -f "$vcfFile" ]; then
		echo "File not found!"
		continue
	fi
	
	#2. Get the cns.call file
	cnFile=`ls "$d"/output/*.cns.call`
	if [ ! -f "$cnFile" ]; then
		echo "File not found!"
		continue
	fi
	
	#3. First generate the input data for PyClone from the snv vcf files
	#Get the sample name for the input file such that PyClone uses this sample name in the data
	filename=$(basename -- "$d")
	sampleName="${filename%.*}"
	
	if [ "$sampleName" = "PBL36" ] || [ "$sampleName" = "NAP38" ] || [ "$sampleName" = "PBL39" ] || [ "$sampleName" = "PBL37" ]; then #would be neater if all of these are provided by theuser
		echo "Excluding reference sample $sampleName"
		continue
	fi
		
	outFile="$d"/"$sampleName".tsv
	python generatePyCloneInput.py "$vcfFile" "$cnFile" "$outFile"
	
	#4. Genreate input string with file locations for PyClone
	if [ "$fileNum" -lt 1 ]; then
		pyCloneInputFiles="$outFile"
	else
		pyCloneInputFiles="$pyCloneInputFiles $outFile"
	fi
	fileNum=$((fileNum+1))
done

#4. Run PyClone with the inputs
echo "$pyCloneInputFiles"
echo "$folder/pycloneResults"
#Why is running this command from here not working??? I for now echo it then I can run itn the terminal myself
echo "Run: PyClone run_analysis_pipeline --in_files "$pyCloneInputFiles" --working_dir "$folder/pycloneResults""
#PyClone run_analysis_pipeline --in_files "$pyCloneInputFiles" --working_dir "$folder/pycloneResults"

#3. When all samples have been run with PyClone, parse the data to LICHeE input format
echo "Run: python parsePyCloneOutputToLICHeEInput.py $folder/pycloneResults/tables/loci.tsv $folder/licheeInput.txt"
#python parsePyCloneOutputToLICHeEInput.py "$folder/pycloneResults/tables/loci.tsv" "$folder/licheeInput.txt"

#4. Run LICHeE
#cd ../../../../Methods/lichee/lichee/LICHeE/release/
echo "Run: ./lichee -build -i $folder/licheeInput.txt -o $folder/lichee_output.txt -cp -n 0 -minVAFPresent 0.0001 -maxVAFAbsent 0.9999 -showTree 0 -color -dot -v"
#./lichee -build -i "$folder/licheeInput.txt" -o "$folder/lichee_output.txt" -cp -n 0 -minVAFPresent 0.0001 -maxVAFAbsent 0.9999 -showTree 0 -color -dot -v
