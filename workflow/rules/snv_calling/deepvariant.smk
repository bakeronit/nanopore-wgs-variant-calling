DEEPVARIANT_GPU_sif = config['deepvariant']['gpu']
DEEPVARIANT_CPU_sif = config['deepvariant']['cpu']

## deepvariant only available for R10
rule call_germline_snv_deepvariant:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai="analysis/bam/{flowcell}/{mode}/{sample}.bam.bai"
    output:
        vcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.vcf.gz",
        gvcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.g.vcf.gz"
    log:
        "logs/deepvariant/{flowcell}.{mode}.{sample}.log"
    benchmark:
        "benchmarks/deepvariant/{flowcell}.{mode}.{sample}.benchmark.txt"
    threads: 12 ## use 12 so a 4-GPUs nodes(with 52 CPUs) can run 4 jobs a time, I don't like 13 so.
    resources:
        nvidia_gpu = 1,
        mem = 64,
        walltime = 48
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        singularity run --nv {DEEPVARIANT_GPU_sif} \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type ONT_R104 \
            --ref {input.genome} \
            --reads {input.bam} \
            --sample_name {wildcards.sample} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {threads} | tee -a {log}
        """

rule margin_phasing:
    input:
        bam="analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai="analysis/bam/{flowcell}/{mode}/{sample}.bam.bai",
        genome=config['reference']['file'],
        vcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.vcf.gz"
    output:
        bam="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.haplotagged.bam",
        vcf="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}.phased.vcf"
    params:
        json=lambda w: 'allParams.haplotag.ont-r94g507.json' if w.flowcell == 'R9' else 'allParams.haplotag.ont-r104q20.json',
        prefix="analysis/snvs/deepvariant/{flowcell}/{mode}/{sample}/{sample}"
    log:
        "logs/deepvariant/{flowcell}.{mode}.{sample}.margin.log"
    benchmark:
        "benchmarks/deepvariant/{flowcell}.{mode}.{sample}.margin.benchmark.txt"
    threads: 24
    resources:
        mem = 48,
        walltime=48
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        singularity run {DEEPVARIANT_CPU_sif} \
        margin phase \
        {input.bam} {input.genome} {input.vcf} \
        {params.json} -t {threads} -o {params.prefix} | tee -a {log}
        """
