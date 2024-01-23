modkit = config['modkit']

rule bamTobedmethyl:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai",
        genome = config['reference']['file'],
    output:
        gz="analysis/mod/{flowcell}/{mode}/{sample}.bed.gz",
        tbi="analysis/mod/{flowcell}/{mode}/{sample}.bed.gz.tbi"
    params:
        bed = "analysis/mod/{flowcell}/{mode}/{sample}.bed"
    threads: 12
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{flowcell}.{mode}.{sample}.pileup.log"
    benchmark:
        "benchmarks/modkit/{flowcell}.{mode}.{sample}.benchmark.txt"
    envmodules:
        "htslib/1.17"
    shell:
        """
        {modkit} pileup {input.bam} {params.bed} --ref {input.genome} \
        --threads {threads} --combine-strands --cpg --log-filepath {log}
        bgzip {params.bed}

        tabix -p bed {output.gz}
        """

rule modsummary:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai = "analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    output:
        "analysis/mod/{flowcell}/{mode}/{sample}.mod_summary.txt"
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    log:
        "logs/modkit/{flowcell}.{mode}.{sample}.summary.log"
    shell:
        """
        {modkit} summary -t {threads} {input.bam} --log-filepath {log} > {output}
        """

rule bamtohapmod:
    input:
        hp_tagged_bam = "analysis/snvs/clair3/{flowcell}/{mode}/{sample}.haplotagged.bam",
        hp_tagged_bai = "analysis/snvs/clair3/{flowcell}/{mode}/{sample}.haplotagged.bam.bai",
        genome = config['reference']['file'],
    output:
        "analysis/mod/{flowcell}/{mode}/{sample}_1.bed.gz",
        "analysis/mod/{flowcell}/{mode}/{sample}_2.bed.gz",
        "analysis/mod/{flowcell}/{mode}/{sample}_1.bed.gz.tbi",
        "analysis/mod/{flowcell}/{mode}/{sample}_2.bed.gz.tbi",
        "analysis/mod/{flowcell}/{mode}/{sample}_ungrouped.bed",
    params:
        outdir = "analysis/mod/{flowcell}/{mode}/",
        hp1_bed = "analysis/mod/{flowcell}/{mode}/{sample}_1.bed",
        hp2_bed = "analysis/mod/{flowcell}/{mode}/{sample}_2.bed",
    threads: 8
    resources:
        mem = 10,
        walltime = 12
    envmodules:
        "htslib/1.17"
    log:
        "logs/modkit/{flowcell}.{mode}.{sample}.pileup.HP.log"
    shell:
        """
        {modkit} pileup {input.hp_tagged_bam} {params.outdir} --ref {input.genome} --partition-tag HP \
        --threads {threads} --combine-strands --cpg --prefix {wildcards.sample} --log-filepath {log}
        bgzip {params.hp1_bed}
        bgzip {params.hp2_bed}

        tabix -p bed {output[0]}
        tabix -p bed {output[1]}
        """