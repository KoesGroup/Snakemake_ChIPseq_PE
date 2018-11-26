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
        "../envs/trimmomatic.yaml"
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
        "../envs/fastqc.yaml"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} &>{log}"

rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"Indexing Reference genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    conda:
        "../envs/samtools_bowtie.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        mapped          = WORKING_DIR + "mapped/{sample}.bam",
        unmapped        = [WORKING_DIR + "unmapped/{sample}.fq." + str(i) +".gz" for i in range(1,2)]
    message: "Mapping files {wildcards.sample}"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome",
        unmapped        = WORKING_DIR + "unmapped/{sample}.fq.gz"
    threads: 10
    conda:
        "../envs/samtools_bowtie.yaml"
    log:
        RESULT_DIR + "logs/bowtie/{sample}.log"
    shell:
        """
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.forward} -2 {input.reverse} -U {input.forwardUnpaired},{input.reverseUnpaired} --un-conc-gz {params.unmapped} | samtools view -Sb - > {output.mapped} 2>{log} 
        """    

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        temp(RESULT_DIR + "mapped/{sample}.sorted.bam")
    message:"Sorting {wildcards.sample} bam file"
    threads: 10
    log:
        RESULT_DIR + "logs/samtools/{sample}.sort.log"
    conda:
        "../envs/samtools.yaml"
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
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools rmdup {input} {output.bam} &>{log}
        samtools index {output.bam}
        """
        #samtools manual says "This command is obsolete. Use markdup instead

# rule bedgraph:
#     input:
#         RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
#     output:
#         RESULT_DIR + "bedgraph/{sample}.sorted.rmdup.bedgraph"
#     params:
#         genome = WORKING_DIR + "genome"
#     message:
#         "Creation of {wildcards.sample} bedgraph file"
#     log:
#         RESULT_DIR + "logs/deeptools/{sample}.sorted.rmdup.bedgraph.log"
#     conda:
#         "../envs/bedtools_env.yaml"
#     shell:
#         "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"
