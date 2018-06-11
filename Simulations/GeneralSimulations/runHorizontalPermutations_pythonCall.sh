#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=20G
#$ -l h_rt=2:00:00
#$ -e noise0.03Mu90_err
#$ -o noise0.03Mu90_out

#Run the simulations again on the previous simulations, but then do a horizontal shuffling of the LAF measurements. 

mu="$1" #select 11% because this equals a mu of 10%, as the task ID starts from 0 (so the tumor fraction is 90%)
uuid="$2"
#Loop through the folders and provide the existing UUIDs

python runSimulation.py "$mu" "$uuid"

