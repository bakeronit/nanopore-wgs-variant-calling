script="/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/snakemake-workflow/scripts"
rule bc_summary:
    input:
        "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam"
    output:
        summary = "analysis/qc/basecalling/{flowcell}/{mode}/{sample}/{run}.summary.txt",
        dy = "analysis/qc/basecalling/{flowcell}/{mode}/{sample}/{run}.data_yield.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 4
    shell:
        """
        {DORADO} summary {input} > {output.summary}
        {script}/data_yield.awk {output.summary} > {output.dy}
        """

rule bam_stats:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.mosdepth.global.dist.txt",
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.mosdepth.summary.txt",
    params:
        prefix = "analysis/qc/bam/{flowcell}/{mode}/{sample}"
    threads: 4
    resources:
        mem = 10,
        walltime = 4
    envmodules:
        "mosdepth/0.2.9"
    shell:
        """
        mosdepth -n -t {threads} {params.prefix} {input}
        """

rule bam_cov:
    input:
        "analysis/bam/{flowcell}/{mode}/{sample}.bam"
    output:
        "analysis/qc/bam/{flowcell}/{mode}/{sample}.bamcov.txt"
    threads: 1
    resources:
        mem = 10,
        walltime = 4
    shell:
        """
        /mnt/backedup/home/jiaZ/working/local/bamcov/bamcov {input} -o {output}
        """
