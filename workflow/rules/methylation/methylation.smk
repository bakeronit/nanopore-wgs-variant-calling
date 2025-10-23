rule bamTobedmethyl:
    input:
        bam = "analysis/bam/{sample}.bam",
        bai = "analysis/bam/{sample}.bam.bai",
        genome = config['reference']['file'],
    output:
        gz="analysis/mod/{sample}.bed.gz",
        tbi="analysis/mod/{sample}.bed.gz.tbi"
    params:
        bed = "analysis/mod/{sample}.bed"
    threads: 12
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{sample}.pileup.log"
    benchmark:
        "benchmarks/modkit/{sample}.benchmark.txt"
    shell:
        """
        modkit pileup {input.bam} {params.bed} --ref {input.genome} \
        --threads {threads} --combine-strands --cpg --log-filepath {log}
        bgzip {params.bed}

        tabix -p bed {output.gz}
        """

rule modsummary:
    input:
        bam = "analysis/bam/{sample}.bam",
        bai = "analysis/bam/{sample}.bam.bai"
    output:
        "analysis/mod/{sample}.mod_summary.txt"
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{sample}.summary.log"
    shell:
        """
        modkit summary -t {threads} {input.bam} --log-filepath {log} > {output}
        """

rule bamtohapmod:
    input:
        hp_tagged_bam = "analysis/bam/{sample}.haplotagged.bam",
        hp_tagged_bai = "analysis/bam/{sample}.haplotagged.bam.bai",
        genome = config['reference']['file'],
    output:
        "analysis/mod/{sample}_1.bed.gz",
        "analysis/mod/{sample}_2.bed.gz",
        "analysis/mod/{sample}_1.bed.gz.tbi",
        "analysis/mod/{sample}_2.bed.gz.tbi",
        "analysis/mod/{sample}_ungrouped.bed",
    params:
        outdir =  "analysis/mod/",
        hp1_bed = "analysis/mod/{sample}_1.bed",
        hp2_bed = "analysis/mod/{sample}_2.bed",
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{sample}.pileup.HP.log"
    shell:
        """
        modkit pileup {input.hp_tagged_bam} {params.outdir} --ref {input.genome} --partition-tag HP \
        --threads {threads} --combine-strands --cpg --prefix {wildcards.sample} --log-filepath {log}
        bgzip {params.hp1_bed}
        bgzip {params.hp2_bed}

        tabix -p bed {output[0]}
        tabix -p bed {output[1]}
        """