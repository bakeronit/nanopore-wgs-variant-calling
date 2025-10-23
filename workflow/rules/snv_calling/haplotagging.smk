Clair3_sif = config['clair3']['sif']

rule haplotagging_whatshap:
    input:
        genome = config['reference']['file'],
        vcf = get_phased_vcf,
        tbi = lambda wildcards: get_phased_vcf(wildcards) + ".tbi",
        bam = "analysis/bam/{sample}.bam"
    output:
        bam = "analysis/bam/{sample}.haplotagged.bam"
    log:
        "logs/whatshap/{sample}.haplotagging_whatshap.log"
    benchmark:
        "benchmarks/whatshap/{sample}.haplotagging_whatshap.benchmark.txt"
    threads: 8
    container: Clair3_sif
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        whatshap haplotag -o {output.bam} \
        --reference {input.genome} {input.vcf} {input.bam} \
        --ignore-read-groups --tag-supplementary --skip-missing-contigs \
        --output-threads={threads} &> {log}
        """

rule index_haplotagged_bam:
    input:
        "analysis/bam/{sample}.haplotagged.bam",
    output:
        "analysis/bam/{sample}.haplotagged.bam.bai",
    threads: 8
    resources:
        mem = 10,
        walltime = 2
    shell:
        """
        samtools index -@8 {input}
        """