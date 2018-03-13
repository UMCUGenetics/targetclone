#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=12G
#$ -l h_rt=10:00:00
#$ -e permutation_err
#$ -o permutation_out

uuid="$1"


python runPermutation.py "$uuid"
