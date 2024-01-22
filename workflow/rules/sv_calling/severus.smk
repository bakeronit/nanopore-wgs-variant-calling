### somatic SV caller for long reads. works better with phased vcf and bam
## it can do multi-sample somatic SV calling, but this workflow only works for single sample (with matched normal) somatic SV calling.
## input: phased normal/tumour bam files, and phased vcf for normal sample.
## need to specify which phased result should be used here, either clair3, pepper, deepvariant.

rule call_somatic_sv_severus:
    input:
        hp_tagged_tumour_bam = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_t}.haplotagged.bam", 
        hp_tagged_normal_bam = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_n}.haplotagged.bam",
        hp_tagged_tumour_bai = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_t}.haplotagged.bam.bai", 
        phased_normal_bai = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_n}.haplotagged.bam.bai",
        phased_vcf = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_n}/phased_merge_output.vcf.gz",
        vntr_bed = config['severus']['vntr']
    output:
        "analysis/svs/severus/{flowcell}/{mode}/{sample_t}.{sample_n}/somatic_SVs/severus_somatic_{sample_t}.haplotagged.vcf"
    log:
        "logs/severus/{flowcell}.{sample_t}.{sample_n}.{flowcell}.{mode}.log"
    benchmark:
        "benchmarks/severus/{flowcell}.{sample_t}.{sample_n}.{flowcell}.{mode}.txt"
    params:
        outdir = "analysis/svs/severus/{flowcell}/{mode}/{sample_t}.{sample_n}"    
    threads: 20
    resources:
        mem = 64,
        walltime = 48
    envmodules:
        "conda-envs/severus-0.1.1"
    shell:
        """
        severus --target-bam {input.hp_tagged_tumour_bam} \
            --control-bam {input.hp_tagged_normal_bam} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --vntr-bed {input.vntr_bed} \
            --phasing-vcf {input.phased_vcf} &>{log}
        """


rule index_haplotagged_bam:
    input:
        "analysis/snvs/clair3/{flowcell}/{mode}/{sample}.haplotagged.bam",
    output:
        "analysis/snvs/clair3/{flowcell}/{mode}/{sample}.haplotagged.bam.bai"
    threads: 8
    resources:
        mem = 10,
        walltime = 8
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools index -@8 {input}
        """
