#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=20G
#$ -l h_rt=2:00:00
#$ -e noise0.03Mu90
#$ -o noise0.03Mu90

#Run a simulation with horizontal dependency removed permutation. We do this 100 times for 1 mu and 1 noise level. 

mu=11 #select 11% because this equals a mu of 10%, as the task ID starts from 0 (so the tumor fraction is 90%)

uuid=$(uuidgen)

python runSimulation.py "$mu" "$uuid"
