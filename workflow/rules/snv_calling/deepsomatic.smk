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
        model="/mnt/backedup/home/jiaZ/working/data/ont_models/dpsomatic_model/weights-143-0.987994.ckpt"
    log:
        "logs/deepsomatic/{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/deepsomatic/{sample_t}.{sample_n}.benchmark.txt"
    threads: 24
    envmodules:
        "singularity/3.7.1"
    resources:
        mem=30,
        walltime=200
    shell:
        """
        singularity exec {DP_somatic_sif} /opt/deepvariant/bin/deepsomatic/run_deepsomatic  \
        --model_type=ONT \
        --ref={input.genome} \
        --reads_normal={input.normal_bam} \
        --reads_tumor={input.tumor_bam} \
        --output_vcf={output} \
        --sample_name_tumor="{wildcards.sample_t}" \
        --sample_name_normal="{wildcards.sample_n}" \
        --num_shards={threads} \
        --logging_dir={log} 
        """
