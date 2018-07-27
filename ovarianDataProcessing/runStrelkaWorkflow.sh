#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=60:00:00
#$ -e strelka_workflow_e
#$ -o strelka_workflow__o


workFlowPath="$1"
jobNum="$2"

python "$workFlowPath"/runWorkflow.py -m sge -j "$jobNum" 

