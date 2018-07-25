# Snakemake file for ChIP-Seq analysis

###############
# Libraries
###############
import os

###############
# Configuration
###############
configfile: "config_tomato.yaml"

WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

GENOME_FASTA_URL = config["refs"]["genome_url"]
GENOME_FASTA_FILE = os.path.basename(config["refs"]["genome_url"])


################
# Desired output
################
rule all:
	input:
		  expand("mapped/{sample}.bam", sample=config["samples"])
	message: "ChIP-seq pipeline succesfully run. The {WORKING_DIR} working directory was removed"
	shell:"rm -rf {WORKING_DIR}"


#############
# Functions #
#############


###############
# Rules
###############
rule get_genome_fasta:
    output:
        WORKING_DIR + "genome.fasta"
    message:"downloading {GENOME_FASTA_FILE} genomic fasta file"
    shell: "wget -O {output} {GENOME_FASTA_URL}"

rule trimmomatic:
    input:
        forward_reads = lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reverse_reads = lambda wildcards: config["samples"][wildcards.sample]["reverse"],
        adapters = config["adapters"]
    output:
        forward = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}.log"
    params :
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward} {output.forwardUnpaired} "
        "{output.reverse} {output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule bowtie2:
    input:
	    forward = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
		forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz"
    output:
	    "mapped/{sample}.bam"
	message: "Mapping files"
	params:
	    bowtie = " ".join(config["bowtie2"]["params"].values())

    threads: 10
    shell:
        "bowtie2 --end-to-end --very-sensitive -X 600 -I 50 "
        "-p {threads} -q --mm -x ../../bowtie/hg19 "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} |"
        "samtools view -Sb {output}"
	# 	check the use of -X and -I for munimum and maximum fragment length for valid paired-end alignments.
	#	-X 500 maximum fragment length for valid paired-end alignments
	#	-I 80 The minimum fragment length for valid paired-end alignments


