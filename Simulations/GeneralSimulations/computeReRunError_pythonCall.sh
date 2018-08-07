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

folder="$1"
prefix="$2"

python visualizeRerunError.py "$folder" "$prefix"
