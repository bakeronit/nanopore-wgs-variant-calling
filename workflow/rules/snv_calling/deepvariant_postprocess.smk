rule call_germline_snv_deepvariant_vcf_stats_report:
    input:
        "analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz",
    output:
        "analysis/snvs/deepvariant/{sample}/{sample}.visual_report.html"
    params:
        outbase = "analysis/snvs/deepvariant/{sample}/{sample}"
    threads: 1
    resources:
        mem = 10,
        walltime = 10
    container: DEEPVARIANT_CPU_sif
    shell:
        """
        vcf_stats_report \
        --input_vcf {input} \
        --outfile_base {params.outbase}
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