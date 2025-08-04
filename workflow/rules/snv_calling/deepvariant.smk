DEEPVARIANT_CPU_sif = config['deepvariant']['cpu']
PEPPER_MARGIN_DV_sif = config['pepper']['sif']

rule call_germline_snv_deepvariant:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{sample}.bam",
        bai="analysis/bam/{sample}.bam.bai"
    output:
        vcf="analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz",
        gvcf="analysis/snvs/deepvariant/{sample}/{sample}.g.vcf.gz"
    log:
        "logs/deepvariant/{sample}.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.benchmark.txt"
    threads: 30 ## use 12 so a 4-GPUs nodes(with 52 CPUs) can run 4 jobs a time.
    resources:
        mem = 64,
        walltime = 48
    container: DEEPVARIANT_CPU_sif
    shell:
        """
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
        bam="analysis/bam/{sample}.bam",
        bai="analysis/bam/{sample}.bam.bai",
        genome=config['reference']['file'],
        vcf="analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz"
    output:
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.phased.vcf.gz"
    params:
        prefix="analysis/snvs/deepvariant/{sample}/{sample}",
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.phased.vcf",
        json="/opt/margin_dir/params/phase2/allParams.haplotag.ont-r104q20.json" #haplotag mode.
    log:
        "logs/deepvariant/{sample}.margin_phasing.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.margin_phasing.benchmark.txt"
    threads: 24
    resources:
        mem = 48,
        walltime=48
    shell:
        """
        margin phase \
        {input.bam} {input.genome} {input.vcf} \
        {params.json} -t {threads} -o {params.prefix} --skipHaplotypeBAM | tee -a {log}

        bgzip {params.vcf}
        tabix -p vcf {output.vcf}
        """
