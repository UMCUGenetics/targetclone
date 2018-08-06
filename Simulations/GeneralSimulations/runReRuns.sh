#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=20G
#$ -l h_rt=2:00:00
#$ -e reruns_e
#$ -o reruns_o

#Run the simulations again on the previous simulations, but then do a horizontal shuffling of the LAF measurements. 

mu=11 #select 11% because this equals a mu of 10%, as the task ID starts from 0 (so the tumor fraction is 90%)
folder="$1"

#Go through the folder. For each subFolder, we need to re-run this 100 times.

#For each folder, start an array job which will run TC on the same dataset 100 times.
#Wait before submitting more jobs.

counter=1
for d in "$folder"/*/ ; do
    
	uuid="$(basename "$d")"
	qsub -t 1-101:1 runReRuns_pythonCall.sh "$mu" "$uuid"
	
	#Per 300 jobs, wait for everything of the previous iteration to complete.
	
	if ! ((counter % 3)); then
		sleep 20m #wait for most of the jobs to complete before re-submitting
	fi
	
	counter=$((counter+1))
done


