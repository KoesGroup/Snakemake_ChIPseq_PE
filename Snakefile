# Snakemake file for ChIP-Seq PE analysis
# The pipeline runs with the command "snakemake --use-conda [--core] ", within the folder containing the Snakefile, after activation of the global_env.yaml

###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version
#############################################
# Configuration and sample sheets
#############################################

configfile: "config.yaml"

WORKING_DIR         = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL    = config["refs"]["genome_url"]
GENOME_FASTA_FILE   = os.path.basename(config["refs"]["genome_url"])
GFF_URL             = config["refs"]["gff_url"]
GFF_FILE            = os.path.basename(config["refs"]["gff_url"])

###########
# Container
###########
singularity: "shub://truatpasteurdotfr/singularity-docker-miniconda"

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_samples_per_treatment(input_df="units.tsv",colsamples="sample",coltreatment="condition",treatment="control"):
    """This function returns a list of samples that correspond to the same experimental condition"""
    df = pd.read_table(input_df)
    df = df.loc[df[coltreatment] == treatment]
    filtered_samples = df[colsamples].tolist()
    return filtered_samples

########################
# Samples and conditions
########################

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

CASES = get_samples_per_treatment(treatment="treatment")
CONTROLS = get_samples_per_treatment(treatment="control")

GROUPS = {
    "group1" : ["ChIP1", "ChIP2", "ChIP3"],
    "group2" : ["ChIP4", "ChIP5", "ChIP6"]
}                                           #I used this dictionnary to define the group of sample used in the multiBamSummary, might be improved a lot

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##################
# Desired output
##################

FASTQC_REPORTS  =     expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"})
BIGWIG          =     expand(RESULT_DIR + "bigwig/{sample}.bw", sample=SAMPLES)
BED_NARROW      =     expand(RESULT_DIR + "bed/{sample}_peaks.narrowPeak", sample=SAMPLES)
MULTIBAMSUMMARY =     RESULT_DIR + "multiBamSummary/MATRIX.npz"
PLOTCORRELATION =     RESULT_DIR + "plotCorrelation/MATRIX.png"
COMPUTEMATRIX   =     expand(RESULT_DIR + "computematrix/{sample}.{type}.gz", sample=SAMPLES, type={"TSS", "scale-regions"})
HEATMAP         =     expand(RESULT_DIR + "heatmap/{sample}.{type}.pdf", sample=SAMPLES, type={"TSS", "scale-regions"})
PLOTFINGERPRINT =     RESULT_DIR + "plotFingerprint/Fingerplot.pdf"
PLOTPROFILE_PDF =     expand(RESULT_DIR + "plotProfile/{sample}.{type}.pdf", sample=SAMPLES, type={"TSS", "scale-regions"})
PLOTPROFILE_BED =     expand(RESULT_DIR + "plotProfile/{sample}.{type}.bed", sample=SAMPLES, type={"TSS", "scale-regions"})
MULTIQC         =     "qc/multiqc.html"

###############
# Final output
################
rule all:
    input:
        FASTQC_REPORTS,
        BIGWIG,
        BED_NARROW,
        MULTIBAMSUMMARY,
        PLOTCORRELATION,
        COMPUTEMATRIX,
        HEATMAP,
        PLOTFINGERPRINT,
        PLOTPROFILE_PDF,
        PLOTPROFILE_BED,
        #MULTIQC
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"rm -rf {WORKING_DIR}"

###############
# Rules
###############

include : "rules/external_data.smk"
include : 'rules/pre_processing.smk'
include : "rules/macs2_peak_calling.smk"
include : "rules/deeptools_post_processing.smk"


rule multiqc:
    input:
        expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"}),
        expand(RESULT_DIR + "bed/{sample}_peaks.xls", sample= SAMPLES)
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"
