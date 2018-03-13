#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=12G
#$ -l h_rt=2:00:00
#$ -e noise2Mu10Clones2
#$ -o noise2Mu10Clones2

#Run a simulation with horizontal dependency removed permutation. We do this 100 times for 1 mu and 1 noise level. 

mu=11 #select 11% because this equals a mu of 10%, as the task ID starts from 0

uuid=$(uuidgen)

python runSimulation.py "$mu" "$uuid"

#qsub runPermutation.sh "$uuid" #disable running permutations for now
