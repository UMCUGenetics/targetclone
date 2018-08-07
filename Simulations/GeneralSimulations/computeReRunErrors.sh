#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=2G
#$ -l h_rt=2:00:00
#$ -e rerunsErrorComputation_e
#$ -o rerunsErrorComputation_o

folder="$1"

#search through the folder for all unique prefixes

simulationDirs=`find "$folder" -maxdepth 1 -type d -not -name "*_*"`

find "$folder" -maxdepth 1 -type d -not -name "*_*" | while read d; do
	
	if [ "$d" == "$folder" ]
	then
		continue
	fi

	uuid="$(basename "$d")"
	echo "$uuid"
	
	#Call the error computation script for this simulation set
	qsub computeReRunError_pythonCall.sh "$folder" "$uuid"
	
done