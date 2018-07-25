# DMC1_ChIPseq pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)


# Scripts

## ChipSeq_BamFilter.py
* *ChipSeq_BamFilter.py* is a python script that will adapt the bam filtering pipeline for each bam file 


# Usage 

## Conda environment
First, you need to install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html).  
1. Create a virtual environment named "dmc1" from the `environment.yaml` file with the following command: `conda env create --name dmc1 --file environment.yaml`
2. Then, activate this virtual environment with `source activate dmc1`    
Now, all softwares and packages versions in use are the one listed in the `environment.yaml` file.

## Configuration file
You will need to provide 

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order. 
It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Dry run
Use the command `snakemake -np` to perform a dry run that prints out the rules and commands. 
 
## Real run
Simply type `Snakemake` and provide the number of cores with `--cores 10` for ten cores for instance.  
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution). 

# Main outputs


