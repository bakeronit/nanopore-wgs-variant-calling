DP_somatic_sif = config['deepsomatic']['cpu']

primary_chromosomes = [str(i) for i in range(1,23)] + ['X', 'Y', 'M']
chromosomes = [f"chr{c}" for c in primary_chromosomes]

wildcard_constraints:
    chr = "|".join(chromosomes)

rule call_somatic_snv_deepsomatic_by_chr:
    """
    separate analysis by chr can avoid requesting big memory
    """
    input:
        genome=config['reference']['file'],
        tumor_bam = "analysis/bam/{sample_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.{chr}.vcf.gz"
    params:
        logdir = "logs/deepsomatic/{sample_t}.{sample_n}.{chr}"
    benchmark:
        "benchmarks/deepsomatic/{sample_t}.{sample_n}.{chr}.benchmark.txt"
    threads: 24
    container: DP_somatic_sif
    resources:
        mem=32,  # some samples need more memory
        walltime=12
    shell:
        """
        run_deepsomatic \
        --model_type=ONT \
        --ref={input.genome} \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor={wildcards.sample_t} \
        --sample_name_normal={wildcards.sample_n} \
        --regions={wildcards.chr} \
        --num_shards={threads} \
        --logging_dir={params.logdir}
        """

rule merge_vcfs_by_chr:
    input:
        [f"analysis/snvs/deepsomatic/{{sample_t}}.{{sample_n}}/output.{chr}.vcf.gz" for chr in chromosomes]
    output:
        vcf = "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.somatic.vcf.gz",
        tbi = "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.somatic.vcf.gz.tbi"
    threads: 1
    shell:
        """
        bcftools concat {input} | \
            bcftools view -i 'GT="1/1"' -f PASS | \
            bcftools sort -Oz -o {output.vcf} 

        bcftools index -t {output.vcf}
        """

