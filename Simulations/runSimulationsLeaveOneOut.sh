#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=12G
#$ -l h_rt=02:00:00
#$ -e leaveOutAll_err
#$ -o leaveOutAll_out

#Run 10 iterations for a simulation where we each time leave out another 

leaveOut=1
mu=26 #-1 will give a mu of 75%, this is the normal fraction

qsub -t 1-101:1 -tc 50 runSimulationLeaveOut.sh "$mu" "$leaveOut"
