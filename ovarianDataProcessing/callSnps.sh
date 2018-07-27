#!/bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=60:00:00
#$ -e snpCalling_e
#$ -o snpCalling_o

bamFile="$1"
outFile="$2"

samtools mpileup --skip-indels -f /hpc/cog_bioinf/ridder/users/mnieboer/from_hpc_tmp/mnieboer/ref/hg19.fa -g "$bamFile" | bcftools call --skip-variants indels --variants-only -m -Oz -o "$outFile"

