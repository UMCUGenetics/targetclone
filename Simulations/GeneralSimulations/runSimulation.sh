#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=02:00:00
#$ -e simulations_e
#$ -o simulations_o

mu="$SGE_TASK_ID"
uuid=$(uuidgen)

python runSimulation.py "$mu" "$uuid"
