### somatic SV caller for long reads. works better with phased vcf and bam
## it can do multi-sample somatic SV calling, but this workflow only works for single sample (with matched normal) somatic SV calling.
## the input will be phased first.

rule call_somatic_sv_severus:
    input:
        phased_tumour_bam = 
        phased_normal_bam = 
        phased_vcf = 
        vntr_bed = config['severus']['vntr']
    output:
        "analysis/sv/severus/{flowcell}/{mode}/{sample}/somatic_SVs/severus_somatic_{sample_t}.haplotagged.vcf"
    params:
        outdir = "analysis/sv/severus/{flowcell}/{mode}/{sample}"    
    threads: 20
    envmodules:
        "conda-envs/severus-0.1.1"
    shell:
        """
        severus --target-bam {input.phased_tumour_bam} \
            --control-bam {input.phased_normal_bam} \
            --out-dir {params.dir} \
            --threads {threads} \
            --vntr-bed {input.vntr_bed} \
            --phasing-vcf {input.phased_vcf}
        """