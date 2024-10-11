DEEPVARIANT_GPU_sif = config['deepvariant']['gpu']
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
        #nvidia_gpu = 1,
    envmodules:
        "singularity/3.7.1"
    shell:
        """
        singularity run --nv {DEEPVARIANT_CPU_sif} \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type ONT_R104 \
            --ref {input.genome} \
            --reads {input.bam} \
            --sample_name {wildcards.sample} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {threads} | tee -a {log}
        """

rule call_germline_snv_deepvariant_s1_make_examples:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{sample}.bam",
        bai="analysis/bam/{sample}.bam.bai"
    output:
        examples="analysis/snvs/deepvariant/{sample}/examples/make_examples.tfrecord.gz"
    params:
        examples_path = "analysis/snvs/deepvariant/{sample}/examples"
    log:
        "logs/deepvariant/{sample}.make_examples.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.make_examples.benchmark.txt"
    threads: 30
    resources:
        mem = 64,
        walltime = 48
    envmodules:
        "parallel/20151122",
        "singularity/3.7.1"
    shell:
        """
        seq 0 $(( {threads} - 1 )) | parallel -q --halt 2 --line-buffer \
        singularity run {DEEPVARIANT_CPU_sif} /opt/deepvariant/bin/make_examples \
            --mode calling \
            --ref {input.genome} \
            --reads {input.bam} \
            --examples "{params.examples_path}/make_examples.tfrecord@{threads}.gz" \
            --add_hp_channel --alt_aligned_pileup "diff_channels" \
            --gvcf {output.examples}.gvcf.tfrecord@{threads}.gz \
            --max_reads_per_partition 600  --min_mapping_quality 5 \
            --parse_sam_aux_fields --partition_size 25000 --phase_reads \
            --pileup_image_width 199 --norealign_reads \
            --sample_name {wildcards.sample} --sort_by_haplotypes \
            --track_ref_reads --vsc_min_fraction_indel 0.12 --vsc_min_fraction_snps 0.08 \
            --task {}
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
        json=lambda w: '/opt/margin_dir/params/phase/allParams.haplotag.ont-r94g507.json' if w.flowcell == 'R9' else '/opt/margin_dir/params/phase2/allParams.haplotag.ont-r104q20.json',
        prefix="analysis/snvs/deepvariant/{sample}/{sample}",
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.phased.vcf",
        params_path = config['pepper']['params']
    log:
        "logs/deepvariant/{sample}.margin_phasing.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.margin_phasing.benchmark.txt"
    threads: 24
    resources:
        mem = 48,
        walltime=48
    envmodules:
        "singularity/3.7.1",
        "htslib/1.16"
    shell:
        """
        singularity exec --bind {params.params_path}:/opt/margin_dir/params/phase2 {PEPPER_MARGIN_DV_sif} \
        margin phase \
        {input.bam} {input.genome} {input.vcf} \
        {params.json} -t {threads} -o {params.prefix} --skipHaplotypeBAM | tee -a {log}

        bgzip {params.vcf}
        tabix -p vcf {output.vcf}
        """
