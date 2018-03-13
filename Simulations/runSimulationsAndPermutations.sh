#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=12G
#$ -l h_rt=15:00:00
#$ -e noise0_err
#$ -o noise0_out

#Run 1 simulation and 10 permutations, split across nodes
#We run 1 simulation for each mu, then we bin these values

mu="$SGE_TASK_ID"

uuid=$(uuidgen)

python runSimulation.py "$mu" "$uuid"

for i in {1..10}
do
   qsub runPermutation.sh "$uuid"
done