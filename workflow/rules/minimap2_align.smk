# TODO: make NanoFlit optional in the future.
rule align_minimap2:
    input:
        ubam = "analysis/ubam/{flowcell}/{mode}/{sample}/{run}.ubam",
        genome = config['reference']['file']
    output:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}/{run}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}/{run}.bam.bai"
    envmodules:
        "samtools/1.17",
        "minimap2/2.26",
        "python/3.9.13"   ## NanoFilt install within this module.
    threads: 24
    resources:
        mem = 36,
        walltime = 2
    params:
        qs = config['filter']['read_qs'],
        headcrop = config['filter']['headcrop'],
        tailcrop = config['filter']['tailcrop'],
        minlen = config['filter']['minlen'],
        rg = "\"@RG\\tID:{sample}.{run}\\tPL:ONT\\tSM:{sample}\""
    benchmark:
        "benchmarks/minimap2/{flowcell}.{sample}.{run}.{mode}.benchmark.txt"
    log:
        "logs/minimap2/{flowcell}.{sample}.{run}.{mode}.log"
    shell: 
        """
        samtools view -e '[qs] >= {params.qs}' {input.ubam} | \
            samtools fastq -@8 -T "*" | \
            NanoFilt --headcrop {params.headcrop} --tailcrop {params.tailcrop} -l {params.minlen} | \
            minimap2 -R {params.rg} \
            -y --MD -ax map-ont -t {threads} {input.genome} - | \
            samtools sort -@8 -O BAM --write-index -o {output.bam}##idx##{output.bai} - &>{log}
        """ 

rule merge_bam:
    input:
        get_bam_of_runs
    output:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    threads: 10
    resources:
        mem = 20,
        walltime = 1
    run:
        if len(input)>1:
            shell("module load samtools/1.17 && samtools merge -@{threads} --write-index {output.bam}##idx##{output.bai} {input}")
        else:
            shell("ln -sr {input} {output.bam} && ln -sr {input}.bai {output.bai} && touch -h {output.bam} && touch -h {output.bai}")
