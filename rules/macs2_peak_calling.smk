rule call_narrow_peaks:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{sample}_peaks.narrowPeak"
    message:
        "Calling narrowPeak for {wildcards.sample}"
    log:RESULT_DIR + "logs/macs2/{sample}_peaks.narrowPeak.log"
    params:
        name        = "{sample}",
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue']),
        outdir      = str(config['macs2']['outdir'])
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir {params.outdir} &>{log}
        """
