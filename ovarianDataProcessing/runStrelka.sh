#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=1:00:00
#$ -e strelka_config_e
#$ -o strelka_config__o



strelkaPath="$1"
normalBam="$2"
tumorBam="$3"
ref="$4"
outDir="$5"
callRegions="$6"

"$strelkaPath"/bin/configureStrelkaSomaticWorkflow.py \
--normalBam "$normalBam" \
--tumorBam "$tumorBam" \
--referenceFasta "$ref" \
--runDir "$outDir" \
--callRegions "$callRegions"
