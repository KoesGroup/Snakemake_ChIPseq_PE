# Snakemake file for ChIP-Seq PE analysis

###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version
#############################################
# Configuration and sample sheets
#############################################

configfile: "configs/config_tomato_sub.yaml"

WORKING_DIR         = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL    = config["refs"]["genome_url"]
GENOME_FASTA_FILE   = os.path.basename(config["refs"]["genome_url"])
GFF_URL             = config["refs"]["gff_url"]
GFF_FILE            = os.path.basename(config["refs"]["gff_url"])

TOTALCORES          = 16                             #check this via 'grep -c processor /proc/cpuinfo'

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

def is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    return pd.isnull(units.loc[(sample), "fq2"])

##############
# Samples and conditions
##############

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

##############
# Desired output
##############

FASTQC_REPORTS  =     expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"})
BAM_INDEX       =     expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam.bai", sample=SAMPLES)
BAM_RMDUP       =     expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam", sample=SAMPLES)
BEDGRAPH        =     expand(RESULT_DIR + "bedgraph/{sample}.sorted.rmdup.bedgraph", sample=SAMPLES)
BIGWIG          =     expand(RESULT_DIR + "bigwig/{sample}.bw", sample=SAMPLES)
BAM_COMPARE     =     expand(RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw", zip, treatment = CASES, control = CONTROLS) #add zip function in the expand to compare respective treatment and control
BED_NARROW      =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.narrowPeak", zip, treatment = CASES, control = CONTROLS)
BED_BROAD       =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.broadPeak", zip, treatment = CASES, control = CONTROLS)
MULTIBAMSUMMARY =     expand(RESULT_DIR + "multiBamSummary/{group}.npz", group = list(GROUPS.keys()))
PLOTCORRELATION =     expand(RESULT_DIR + "plotCorrelation/{sample}.png", sample = list(GROUPS.keys()))
COMPUTEMATRIX   =     expand(RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz", treatment = CASES, control = CONTROLS)
HEATMAP         =     expand(RESULT_DIR + "heatmap/{treatment}_{control}.pdf", treatment = CASES, control = CONTROLS)
PLOTFINGERPRINT =     expand(RESULT_DIR + "plotFingerprint/{treatment}_vs_{control}.pdf", zip, treatment = CASES, control = CONTROLS)

###############
# Final output
################
rule all:
    input:
        BAM_INDEX,
        BAM_RMDUP,
        FASTQC_REPORTS,
        BEDGRAPH,
        BIGWIG,
        BAM_COMPARE,
        BED_NARROW,
        #BED_BROAD
        MULTIBAMSUMMARY,
        PLOTCORRELATION,
        COMPUTEMATRIX,
        HEATMAP,
        PLOTFINGERPRINT
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"#rm -rf {WORKING_DIR}"

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
        reads = get_fastq,
        adapters = config["adapters"]
    output:
        forward_reads   = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse_reads   = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}.log"
    benchmark:
        RESULT_DIR + "benchmarks/{sample}.trimmomatic.benchmark.txt"
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
    conda:
        "envs/trimmomatic_env.yaml"
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} &>{log}"

rule fastqc:
    input:
        fwd = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        rev = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz"
    output:
        fwd = RESULT_DIR + "fastqc/{sample}_forward_fastqc.zip",
        rev = RESULT_DIR + "fastqc/{sample}_reverse_fastqc.zip"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastqc.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "---Quality check of trimmed {wildcards.sample} sample with FASTQC"
    conda:
        "envs/fastqc_env.yaml"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} &>{log}"

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
    benchmark:
        RESULT_DIR + "benchmarks/genome.index.benchmark.txt"
    threads: 10
    conda:
        "envs/samtools_bowtie_env.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        temp(WORKING_DIR + "mapped/{sample}.bam")
    message: "Mapping files {wildcards.sample}"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome"
    threads: 10
    benchmark:
        RESULT_DIR + "benchmarks/{sample}.align.benchmark.txt"
    conda:
        "envs/samtools_bowtie_env.yaml"
    shell:
        "bowtie2 {params.bowtie} "
        "--threads {threads} "
        "-x {params.index} "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} "   # also takes the reads unpaired due to trimming
        "| samtools view -Sb - > {output}"                       # to get the output as a BAM file directly

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        temp(RESULT_DIR + "mapped/{sample}.sorted.bam")
    message:"sorting {wildcards.sample} bam file"
    benchmark:
        RESULT_DIR + "benchmarks/{sample}.sort.benchmark.txt"
    threads: 10
    log:
        RESULT_DIR + "logs/samtools/{sample}.sort.log"
    conda:
        "envs/samtools.yaml"
    shell:"samtools sort -@ {threads} -o {output} {input} &>{log}"

rule rmdup:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.bam"
    output:
        bam = RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam",
        bai = RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam.bai"        #bai files required for the bigwig and bamCompare rules
    message: "Removing duplicate from file {wildcards.sample}"
    log:
        RESULT_DIR + "logs/samtools/{sample}.sorted.rmdup.bam.log"
    benchmark:
        RESULT_DIR + "benchmarks/{sample}.rmdup.benchmark.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools rmdup {input} {output.bam} &>{log}
        samtools index {output.bam}
        """
        #samtools manual says "This command is obsolete. Use markdup instead

rule bedgraph:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bedgraph/{sample}.sorted.rmdup.bedgraph"
    params:
        genome = WORKING_DIR + "genome"
    message:
        "Creation of {wildcards.sample} bedgraph file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}.sorted.rmdup.bedgraph.log"
    benchmark:
        RESULT_DIR + "benchmarks/{sample}.bedgraph.benchmark.txt"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"

rule bamCoverage:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bigwig/{sample}.bw"
    message:
        "Converting {wildcards.sample} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE     = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]), #take argument separated as a list separated with a space
        EXTENDREADS             = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize                 = str(config['bamCoverage']["params"]['binSize']),
        normalizeUsing          = str(config['bamCoverage']["params"]['normalizeUsing']),
        ignoreForNormalization  = str(config['bamCoverage']["params"]['ignoreForNormalization']),
        smoothLength            = str(config['bamCoverage']["params"]['smoothLength'])
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input} \
        -o {output} \
        --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} \
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        &>{log}"

rule bamcompare:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",              #input requires an indexed bam file
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw"
    message:
        "Running bamCompare for {wildcards.treatment} and {wildcards.control}"
    log:
        RESULT_DIR + "logs/deeptools/log2_{treatment}_{control}.bamcompare.bw.log"
    conda:
        "envs/deeptools.yaml"
    params:
        binSize             = str(config['bamcompare']['binSize']),
        normalizeUsing      = str(config['bamcompare']['normalizeUsing']),
        EFFECTIVEGENOMESIZE = str(config["bamcompare"]["EFFECTIVEGENOMESIZE"]),
        operation           = str(config['bamcompare']['operation']),
        smoothLength        = str(config['bamcompare']['smoothLength']),
        ignoreForNormalization = str(config['bamcompare']['ignoreForNormalization']),
        scaleFactorsMethod  = str(config['bamcompare']['scaleFactorsMethod'])
    shell:
        "bamCompare -b1 {input.treatment} \
        -b2 {input.control}  \
        --binSize {params.binSize} \
        -o {output.bigwig} \
        --normalizeUsing {params.normalizeUsing} \
        --operation {params.operation} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        --scaleFactorsMethod {params.scaleFactorsMethod} \
        &>{log}"

rule call_narrow_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.narrowPeak"
    message:
        "Calling narrowPeak for {wildcards.treatment} and {wildcards.control}"
    params:
        name        = "{treatment}_vs_{control}",        #this option will give the output name, has to be similar to the output
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.narrowPeak.log"
    benchmark:
        RESULT_DIR + "benchmarks/{treatment}_vs_{control}.peaknarrow.benchmark.txt"
    conda:
        "envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/ &>{log}
        """

rule call_broad_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.broadPeak"
    message:
        "Calling broadPeak for {wildcards.treatment} and {wildcards.control}"
    params:
        name        = "{treatment}_vs_{control}",
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.broadPeak.log"
    benchmark:
        RESULT_DIR + "benchmarks/{treatment}_vs_{control}.peakbroad.benchmark.txt"
    conda:
        "envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/ &>{log}
        """


rule multiBamSummary:
    input:
        lambda wildcards: expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam", sample = GROUPS[wildcards.group])
    output:
        RESULT_DIR + "multiBamSummary/{group}.npz"
    message:
        "Computing the read coverage for {wildcards.group} "
    threads: 10
    params:
        binSize     = str(config['multiBamSummary']['binSize'])
    log:
        RESULT_DIR + "logs/deeptools/multibamsummary/{group}.log"

    shell:
        "multiBamSummary bins \
        --bamfiles {input} \
        --numberOfProcessors {threads}\
        --binSize {params.binSize} \
        --centerReads \
        --extendReads \
        -o {output} \
        2> {log}"


rule plotCorrelation:
    input:
        RESULT_DIR + "multiBamSummary/{sample}.npz"
    output:
        RESULT_DIR + "plotCorrelation/{sample}.png"
    log:
        RESULT_DIR + "logs/deeptools/plotcorrelation/{sample}.log"
    params:
        corMethod  = str(config['plotCorrelation']['corMethod']),
        whatToPlot = str(config['plotCorrelation']['whatToPlot']),
        color      = str(config['plotCorrelation']['color'])
    conda:
        "envs/deeptools.yaml"
    shell:
        "plotCorrelation \
                    --corData {input} \
                    --corMethod {params.corMethod} \
                    --whatToPlot {params.whatToPlot} \
                    --skipZeros \
                    --colorMap {params.color} \
                    --plotFile {output} \
                    --plotNumbers \
                    2> {log}"

#--plotTitle 'Pearson Correlation of {params.title} coverage' \

rule get_gff:
    output:
        WORKING_DIR + "gene_model.gff"
    message:"downloading {GFF_FILE} genomic fasta file"
    shell: "wget -O {output} {GFF_URL}"

rule gff_to_gtf:
    input:
        WORKING_DIR + "gene_model.gff"
    output:
        WORKING_DIR + "gene_model.gtf"
    shell:
        "python scripts/gff_to_gtf.py {input} {output}"


rule computeMatrix:
    input:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw",
        bed    = WORKING_DIR + "gene_model.gtf"
    output:
        RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz"
    threads: 10
    params:
        binSize = str(config['computeMatrix']['binSize'])
    conda:
        "envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/computematrix/{treatment}_{control}.log"
    shell:
        "computeMatrix \
        reference-point \
        --referencePoint TSS \
        -S {input.bigwig} \
        -R {input.bed} \
        --afterRegionStartLength 3000 \
        --beforeRegionStartLength 3000 \
        --numberOfProcessors {threads} \
        --binSize {params.binSize} \
        -o {output} \
        2> {log}"

rule plotHeatmap:
    input:
        RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz"
    output:
        RESULT_DIR + "heatmap/{treatment}_{control}.pdf"
    params:
        kmeans = str(config['plotHeatmap']['kmeans']),
        color  = str(config['plotHeatmap']['color']),
        plot   = str(config['plotHeatmap']['plot']),
        cluster = "{treatment}_vs_{control}.bed"
    conda:
        "envs/deeptools.yaml"
    shell:
        "plotHeatmap \
        --matrixFile {input} \
        --outFileName {output} \
        --kmeans {params.kmeans} \
        --colorMap {params.color} \
        --legendLocation best \
        --outFileSortedRegions {params.cluster}"

rule plotFingerprint:
    input:
        treatment = expand(RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam", treatment = CASES),
        control   = expand(RESULT_DIR + "mapped/{control}.sorted.rmdup.bam", control = CONTROLS)
    output:
        pdf = RESULT_DIR + "plotFingerprint/{treatment}_vs_{control}.pdf"
    params:
        EXTENDREADS  = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize      = str(config['bamCoverage']["params"]['binSize'])
    conda:
        "envs/deeptools.yaml"
    shell:
        "plotFingerprint \
        -b {input.treatment} {input.control} \
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --plotFile {output}"

rule plotProfile:
    input:
        RESULT_DIR + "multiBamSummary/{sample}.npz"
    output:
        pdf = RESULT_DIR + "plotProfile/{treatment}_{control}.pdf",
        bed = RESULT_DIR + "plotProfile/{treatment}_{control}.bed"
    params:
        kmeans      = str(config['plotHeatmap']['kmeans']),
        startLabel  = str(config['plotHeatmap']['startLabel']),
        endLabel    = str(config['plotHeatmap']['endLabel'])
    conda:
        "envs/deeptools.yaml"
    shell:
        "plotProfile \
        --matrixFile {input} \
        --outFileName {output.pdf} \
        --outFileSortedRegions {output.bed} \
        --kmeans {params.kmeans} \
        --startLabel {params.startLabel} \
        --endLabel {params.endLabel}"

rule multiqc:
    input:
    output:
        "multiqc_report.html"
    conda:
        "envs/multiqc_env.yaml"
    shell:""
        
