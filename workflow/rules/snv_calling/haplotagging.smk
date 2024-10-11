Clair3_sif = config['clair3']['sif']

rule haplotagging_whatshap:
    input:
        genome = config['reference']['file'],
        vcf = get_phased_vcf,
        bam = "analysis/bam/{sample}.bam"
    output:
        bam = "analysis/bam/{sample}.haplotagged.bam"
    log:
        "logs/whatshap/{sample}.haplotagging_whatshap.log"
    benchmark:
        "benchmarks/whatshap/{sample}.haplotagging_whatshap.benchmark.txt"
    threads: 8
    envmodules:
        "singularity/3.7.1"
    resources:
        mem = 48,
        walltime = 48
    shell:
        """
        singularity exec {Clair3_sif} \
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
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools index -@8 {input}
        """