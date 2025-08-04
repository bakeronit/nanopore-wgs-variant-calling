rule qc_bam_stats:
    input:
        "analysis/bam/{sample}.bam"
    output:
        "analysis/qc/bam/{sample}.mosdepth.global.dist.txt",
        "analysis/qc/bam/{sample}.mosdepth.summary.txt",
    params:
        prefix = "analysis/qc/bam/{sample}"
    threads: 4
    resources:
        mem = 10,
        walltime = 12
    shell:
        """
        mosdepth -n -t {threads} {params.prefix} {input}
        """

rule qc_bam_cov:
    input:
        "analysis/bam/{sample}.bam"
    output:
        "analysis/qc/bam/{sample}.bamcov.txt"
    threads: 1
    resources:
        mem = 32,
        walltime = 12
    shell:
        """
        samtools coverage {input} > {output}
        """