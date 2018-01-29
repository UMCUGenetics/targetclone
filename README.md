# TargetClone
TargetClone: A multi-sample approach for reconstructing subclonal evolution of tumors

## TargetClone

TargetClone is a tool for reconstructing the subclonal evolution of tumors from allele frequency and somatic SNV measurements. It also infers the copy numbers, alleles and fraction from multiple samples of the same tumor. Preferably, the input samples are purified, such as microdissections, and are as homogeneous as possible (containing one tumor subclone and may be contaminated with healthy cells). TargetClone works by assuming that tumor subclones are horizontally and vertically dependent. If an event affects one position, it is very likely that the neighboring position is also affected by the same event (horizontal dependency). Furthermore, subclones inherit these events from their parent, meaning that closely related subclones likely have highly similar mutation patterns as opposed to distantly related subclones (vertical dependency). The infinite sites assumption (ISA) is used to reconstruct the evolutionary tree. If it is known beforehand that the ISA will not hold for the data, it is possible to exclude the respective somatic SNVs (see below).

## Installation instructions

To install TargetClone, users are recommended to download (GitHub zip) and unpack at desired destination. The tool requires the following dependencies (tested with listed versions):

Python 2.7.5 <br />
Numpy 1.7.1 <br />
Scipy 0.12.1 <br />
Bokeh 0.12.5 <br />
Zss 1.1.4 <br />

## Running TargetClone, input and output

After unpacking and installing all required dependencies, TargetClone can be run with desired input data using: 
```
python main.py 'input file.txt' 'output directory/'
```
The input file is assumed to be a tab-delimited file containing the chromosome, position, indication of whether the measurement is a somatic SNV, and the AF measurement at SNP positions (heterozygous in the reference sample) across multiple samples. It is recommended that the measurements of multiple replicates of reference samples are averaged into one sample. Examples of input files are provided in the Examples folder (see below). The output directory will be filled with tab-delimited .txt files containing the inferred C and A matrices, a .txt file listing the tumor fraction (mu) and a .txt file with a text representation of the evolutionary trees. Trees.html contains the visual representation of the tree. 

## Examples

In the Examples folder, two examples are provided to test the output of TargetClone. Example 1 is a simple case with 11 AF measurements and 3 SNVs measured in 3 samples. Example 2 is a full dataset that is part of the simulations (sequencing noise level of 0.02).

To run any example (i.e. example 1), navigate to the directory containing main.py and enter:
```
python testTargetClone.py ../Examples/Example1/Example1.txt ../Examples/Example1/Output/ ../Examples/Example1
```
The output files will be written to the ../Examples/Example1/Output/ folder. Example1.txt contains the AF and SNV measurements. The third argument should point to the location containing the real output data, to which the output will be compared to test the error rates.

## Settings

All settings of TargetClone can be updated in settings.py (found alongside main.py in the TargetClone directory), in which the settings are further detailed.

## Scaling beyond targeted sequencing datasets

On average, TargetClone can complete one run with approximately 8 samples and 450 SNP measurements within two hours on one CPU core. If a dataset contains more than 500 measured SNP positions (for example if the data is generated with exome sequencing), the run time of the algorithm may increase. In these cases, it is recommended to perform a pre-segmentation on the SNP measurements into regions with equal AF. Furthermore, it is recommended to cluster somatic SNV measurements into groups that are shared or absent across samples to reduce any influence of noise. 

## When the infinite sites assumption does not hold

If it is known that the ISA does not hold in the respective dataset, the somatic SNVs that violate the ISA can be excluded from analysis with TargetClone. This file requires the specification of the chromosome and position on which the somatic SNVs are present. An example of such a file is provided in the Examples folder (excludedSNVs.txt) and the location can be specified in settings.py. 
