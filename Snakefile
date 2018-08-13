# Snakemake file for ChIP-Seq analysis

###############
# Libraries
###############
import os

###############
# Configuration
###############
configfile: "config_tomato.yaml"

FASTQ_DIR = config["fastq_dir"]        # where to find the fastq files
WORKING_DIR = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL = config["refs"]["genome_url"]
GENOME_FASTA_FILE = os.path.basename(config["refs"]["genome_url"])

##############
# Wildcards
##############
wildcard_constraints:
    sample="[A-Za-z0-9]+"

################
# Final output
################
rule all:
    input:
        expand(WORKING_DIR + "mapped/{sample}.sorted.bam", sample=config["samples"])
    message: "ChIP-seq pipeline succesfully run. The {WORKING_DIR} working directory was removed"
	#finger crossed to see this message!
    shell:"rm -rf {WORKING_DIR}"


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
        #forward = FASTQ_DIR + "{sample}_1.fastq.gz",
		#reverse = FASTQ_DIR + "{sample}_2.fastq.gz",
		#12/08/18 JC : commented out lines 48,49, it makes more sense for me to call sample with the wildcards
        forward_reads = lambda wildcards: FASTQ_DIR + config["samples"][wildcards.sample]["forward"],
        reverse_reads = lambda wildcards: FASTQ_DIR + config["samples"][wildcards.sample]["reverse"],
        adapters = config["adapters"]
    output:
        forward_reads = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse_reads = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
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
        "{input.forward_reads} "
        "{input.reverse_reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"indexing genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        temp(WORKING_DIR + "mapped/{sample}.bam")
    message: "Mapping files"
    params:
        bowtie = " ".join(config["bowtie2"]["params"].values()),
        index = WORKING_DIR + "genome"
    threads: 10
    shell:
        "bowtie2 {params.bowtie} "
        "--threads {threads} "
        "-x {params.index} "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} "   # also takes the reads unpaired due to trimming
        "|samtools view -Sb - > {output}"                       # to get the output as a BAM file directly

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        WORKING_DIR + "mapped/{sample}.sorted.bam"
    message:"sorting {wildcards.sample} bam file"
    threads: 10
    shell:"samtools sort -@ {threads} -o {output} {input}"
