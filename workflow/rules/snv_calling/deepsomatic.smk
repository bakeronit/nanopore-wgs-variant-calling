DP_somatic_sif = config['deepsomatic']['cpu']

rule call_somatic_snv_deepsomatic:
    input:
        genome=config['reference']['file'],
        tumor_bam = "analysis/bam/{sample_t}.bam",
        tumor_bai = "analysis/bam/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{sample_n}.bam",
        normal_bai = "analysis/bam/{sample_n}.bam.bai"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.vcf.gz"
    params:
        logdir = "logs/deepsomatic/{sample_t}.{sample_n}"
    benchmark:
        "benchmarks/deepsomatic/{sample_t}.{sample_n}.benchmark.txt"
    threads: 24
    container: DP_somatic_sif
    resources:
        mem=36,  # using 30gb will cause Cgroup out of memory error
        walltime=80
    shell:
        """
        run_deepsomatic  \
        --model_type=ONT \
        --ref={input.genome} \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --num_shards={threads} \
        --logging_dir={params.logdir}
        """

rule extract_somatic_snv_deepsomatic:
    input:
        vcf = "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.vcf.gz"
    output:
        "analysis/snvs/deepsomatic/{sample_t}.{sample_n}/output.somatic.vcf.gz"
    threads: 1
    shell:
        """
        bcftools view -i 'GT="1/1"' -f PASS {input.vcf} | bgzip -c > {output}
        tabix -p vcf {output}
        """
