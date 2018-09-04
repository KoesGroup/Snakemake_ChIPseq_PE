# ChIP_seq_Snakemake
A snakemake pipeline for the analysis of ChIP-seq data

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

# Aim
Map paired-end Illumina ChIP-seq data.

# Content
- Snakefile containing the targeted output and the rules to generate them from the input files.
- config/ , folder containing the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules.
- Fastq/, folder containing subsetted paired-end fastq files used to test locally the pipeline. Generated using [Seqtk](https://github.com/lh3/seqtk): `seqtk sample -s100 read1.fq 5000 > sub1.fqseqtk sample -s100 read2.fq 5000 > sub2.fq`
- envs/, folder containing the environment needed for the Snakefile to run. Need to make one specifically for MACS2 as MACS2 uses python 2.7 following the information found [here](https://groups.google.com/forum/#!searchin/snakemake/macs%7Csort:relevance/snakemake/60txGSq81zE/NzCUTdJ_AQAJ).


# Usage

## Conda environment
First, you need to install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html).
1. Create a virtual environment named "chipseq" from the `environment.yaml` file with the following command: `conda env create --name chipseq --file ~/envs/global_env.yaml`
2. Then, activate this virtual environment with `source activate dmc1`
Now, all softwares and packages versions in use are the one listed in the `global_env.yaml` file.

## Configuration file
The `~/configs/config_tomato_sub.yaml` file specifies the sample list, the genomic reference fasta file to use, the directories to use, etc. This file is then used to build parameters in the main `Snakefile`.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order.
It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Samples
Samples are listed in the `units.tsv` and `sample.tsv` files. Change files according to the sample you will use.

## MACS2
Peak calling rules using macs2 require the activation of the 'macs2' environment. Please create the environment before running the snakefile `conda env create --name macs2 --file envs/macs2_env.yaml`. The enviromment is activated within the rules.

## Dry run
Use the command `snakemake -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `Snakemake` and provide the number of cores with `--cores 10` for ten cores for instance.
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).

# Main outputs
The main output are for now sorted bam files, qc files, rmdup.sorted.bam files.
