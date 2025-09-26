"""
To run deepvariant in steps, so the resources on HPC clusters can be managed and release when not needed.
"""
DEEPVARIANT_GPU_sif = config['deepvariant']['gpu']
DEEPVARIANT_CPU_sif = config['deepvariant']['cpu']
## TODO: need to run deepvaraint and check pacbio specific parameters.
#NSHARD = int(config['deepvariant']['make_examples_theads'])
NSHARD=12

def tfrecord_suffix(ncpu: int, json=False):
    suffix = [f".tfrecord-{i:05d}-of-{ncpu:05d}.gz" for i in range(ncpu)]
    if json:
        suffix = [f"{s}.example_info.json" for s in suffix]
    return suffix

rule call_germline_snv_deepvariant_s1_make_examples:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{sample}.bam",
        bai="analysis/bam/{sample}.bam.bai"
    output:
        multiext("analysis/snvs/deepvariant/{sample}/examples/make_examples", *tfrecord_suffix(NSHARD)),
        multiext("analysis/snvs/deepvariant/{sample}/examples/make_examples", *tfrecord_suffix(NSHARD, json=True)),
        multiext("analysis/snvs/deepvariant/{sample}/examples/make_examples_call_variant_outputs", *tfrecord_suffix(NSHARD)),
        multiext("analysis/snvs/deepvariant/{sample}/examples/gvcf", *tfrecord_suffix(NSHARD)),
    params:
        examples_path = "analysis/snvs/deepvariant/{sample}/examples",
        model = "/opt/models/ont_r104",
        smallmodel = "/opt/smallmodels/ont_r104"
    log:
        "logs/deepvariant/{sample}.make_examples.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.make_examples.benchmark.txt"
    threads: NSHARD
    resources:
        mem = 20,
        walltime = 24
    container: DEEPVARIANT_CPU_sif
    shell:
        """
        seq 0 $(( {NSHARD} - 1 )) | parallel -q --halt 2 --line-buffer \
        make_examples \
            --mode calling \
            --ref {input.genome} \
            --reads {input.bam} \
            --examples {params.examples_path}/make_examples.tfrecord@{NSHARD}.gz \
            --checkpoint {params.model} \
            --alt_aligned_pileup diff_channels \
            --call_small_model_examples \
            --gvcf {params.examples_path}/gvcf.tfrecord@{NSHARD}.gz \
            --max_reads_per_partition 600 \
            --min_mapping_quality 5 \
            --parse_sam_aux_fields \
            --partition_size 25000 \
            --phase_reads \
            --pileup_image_width 99 \
            --norealign_reads \
            --sample_name {wildcards.sample} \
            --small_model_indel_gq_threshold 17 \
            --small_model_snp_gq_threshold 9 \
            --small_model_vaf_context_window_size 51 \
            --sort_by_haplotypes \
            --track_ref_reads \
            --trained_small_model_path {params.smallmodel} \
            --trim_reads_for_pileup \
            --vsc_min_fraction_indels 0.12 \
            --vsc_min_fraction_snps 0.08 \
            --task {{}}
        """

rule call_germline_snv_deepvariant_s2_call_variants:
    input:
        rules.call_germline_snv_deepvariant_s1_make_examples.output,
    output:
        "analysis/snvs/deepvariant/{sample}/examples/call_variants_output-00000-of-00001.tfrecord.gz",
    params:
        out_name = "analysis/snvs/deepvariant/{sample}/examples/call_variants_output.tfrecord.gz",
        example_name = lambda w: f"analysis/snvs/deepvariant/{w.sample}/examples/make_examples.tfrecord@{NSHARD}.gz",
        model = "/opt/models/ont_r104"
    log:
        "logs/deepvariant/{sample}.call_variants.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.call_variants.benchmark.txt"
    threads: 1
    resources:
        mem = 20,
        walltime = 24
    container: DEEPVARIANT_CPU_sif
    shell:
        """
        call_variants \
        --outfile {params.out_name} \
        --examples {params.example_name} \
        --checkpoint {params.model}
        """

rule call_germline_snv_deepvariant_s3_postprocess_variants:
    input:
        rules.call_germline_snv_deepvariant_s1_make_examples.output,
        rules.call_germline_snv_deepvariant_s2_call_variants.output,
        ref = config['reference']['file'],
    output:
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz",
        vcf_tbi = "analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz.tbi",
        gvcf = "analysis/snvs/deepvariant/{sample}/{sample}.g.vcf.gz",
        gvcf_tbi = "analysis/snvs/deepvariant/{sample}/{sample}.g.vcf.gz.tbi"
    params:
        input_name = "analysis/snvs/deepvariant/{sample}/examples/call_variants_output.tfrecord.gz",
        cvo_name = lambda w: f"analysis/snvs/deepvariant/{w.sample}/examples/make_examples_call_variant_outputs.tfrecord@{NSHARD}.gz",
        gvcf_name = lambda w: f"analysis/snvs/deepvariant/{w.sample}/examples/gvcf.tfrecord@{NSHARD}.gz",
        model = "/opt/models/ont_r104",
    log:
        "logs/deepvariant/{sample}.postprocess_variants.log"
    benchmark:
        "benchmarks/deepvariant/{sample}.postprocess_variants.benchmark.txt"
    threads: 12
    resources:
        mem = 20,
        walltime = 30
    container: DEEPVARIANT_CPU_sif
    shell:
        """
        postprocess_variants \
        --ref {input.ref} \
        --infile {params.input_name} \
        --outfile {output.vcf} \
        --cpus {threads} \
        --small_model_cvo_records {params.cvo_name} \
        --gvcf_outfile {output.gvcf} \
        --nonvariant_site_tfrecord_path {params.gvcf_name} \
        --sample_name {wildcards.sample}
        """
