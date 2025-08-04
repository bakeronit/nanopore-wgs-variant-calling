### somatic SV caller for long reads. works better with phased vcf and bam
## it can do multi-sample somatic SV calling, but this workflow only works for single sample (with matched normal) somatic SV calling.
## input: phased normal/tumour bam files, and phased vcf for normal sample.
## need to specify which phased result should be used here, either clair3, pepper, deepvariant.

rule call_somatic_sv_severus:
    input:
        hp_tagged_tumour_bam = "analysis/bam/{sample_t}.haplotagged.bam", 
        hp_tagged_normal_bam = "analysis/bam/{sample_n}.haplotagged.bam",
        hp_tagged_tumour_bai = "analysis/bam/{sample_t}.haplotagged.bam.bai", 
        hp_tagged_normal_bai = "analysis/bam/{sample_n}.haplotagged.bam.bai",
        phased_vcf = get_phased_vcf,
        vntr_bed = config['annotation']['vntr']
    output:
        "analysis/svs/severus/{sample_t}.{sample_n}/somatic_SVs/severus_somatic.vcf"
    log:
        "logs/severus/{sample_t}.{sample_n}.log"
    benchmark:
        "benchmarks/severus/{sample_t}.{sample_n}.txt"
    params:
        outdir = "analysis/svs/severus/{sample_t}.{sample_n}",
        min_svlen = 50    
    threads: 20
    resources:
        mem = 64,
        walltime = 48
    shell:
        """
        severus --target-bam {input.hp_tagged_tumour_bam} \
            --control-bam {input.hp_tagged_normal_bam} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --vntr-bed {input.vntr_bed} \
            --min-sv-size {params.min_svlen} \
            --single-bp \
            --between-junction-ins \
            --write-collapsed-dup \
            --output-read-ids \
            --output-LOH \
            --phasing-vcf {input.phased_vcf} &>{log}
        """