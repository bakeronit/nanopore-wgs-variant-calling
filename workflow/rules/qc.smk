#script = "../scripts"
bamcov = config['bamcov']

rule bc_summary:
    input:
        "analysis/ubam/{sample}/{run}.ubam"
    output:
        summary = "analysis/qc/basecalling/{sample}/{run}.summary.txt",
        dy = "analysis/qc/basecalling/{sample}/{run}.data_yield.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 10
    shell:
        """
        {DORADO} summary {input} > {output.summary}
        {script}/data_yield.awk {output.summary} > {output.dy}
        """

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
    envmodules:
        config['modules']['mosdepth']
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
        {bamcov} -H {input} -o {output}
        """