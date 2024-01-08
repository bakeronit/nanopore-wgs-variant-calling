### somatic SV caller for long reads. works better with phased vcf and bam
## it can do multi-sample somatic SV calling, but this workflow only works for single sample (with matched normal) somatic SV calling.
## input: phased normal/tumour bam files, and phased vcf for normal sample.
## need to specify which phased result should be used here, either clair3, pepper, deepvariant.

rule call_somatic_sv_severus:
    input:
        phased_tumour_bam = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_t}.haplotagged.bam", 
        phased_normal_bam = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_n}.haplotagged.bam", 
        phased_vcf = "analysis/snvs/clair3/{flowcell}/{mode}/{sample_n}/phased_merge_output.vcf.gz",
        vntr_bed = config['severus']['vntr']
    output:
        "analysis/svs/severus/{flowcell}/{mode}/{sample_t}.{sample_n}/somatic_SVs/severus_somatic_{sample_t}.{sample_n}.haplotagged.vcf"
    params:
        outdir = "analysis/sv/severus/{flowcell}/{mode}/{sample_t}.{sample_n}"    
    threads: 20
    envmodules:
        "conda-envs/severus-0.1.1"
    shell:
        """
        severus --target-bam {input.phased_tumour_bam} \
            --control-bam {input.phased_normal_bam} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --vntr-bed {input.vntr_bed} \
            --phasing-vcf {input.phased_vcf}
        """