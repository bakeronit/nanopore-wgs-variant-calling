## germline variant calling with pepper-margin-deepvariant pipeline
## it will generates phased vcf and haplotagged bam files
## set --phased_output to enable phasing
## set --skip_final_phased_bam will remove the bam files
##### 

## propably deprecated this, since pepper margin deepvariant haven't been updated for a while.

PEPPER_MARGIN_DV_sif = config['pepper']['sif']

rule call_germline_snv_pepper:
    input:
        genome=config['reference']['file'],
        bam="analysis/bam/{sample}.bam",
        bai="analysis/bam/{sample}.bam.bai"
    output:
        vcf = "analysis/snvs/pepper/{sample}/{sample}.phased.vcf.gz",
        #bam = "analysis/snvs/pepper/{flowcell}/{mode}/{sample}/{sample}.haplotagged.bam" if config['phased_snv_from'] == "pepper" else []
    log:
        "logs/pepper/{sample}.log"
    benchmark:
        "benchmarks/pepper/{sample}.benchmark.txt"
    threads: 24
    resources:
        mem = 64,
        walltime = 100
    container: PEPPER_MARGIN_DV_sif
    params:
        output_dir = "analysis/snvs/pepper/{sample}",
        mode = '--ont_r10_q20', ## TODO: need to be configurable to hifi or R9 ONT, but maybe this pepper_margin_deepvariant will be deprecated.
        #keep = "--skip_final_phased_bam" if config['phased_snv_from'] != "pepper" else "" ## save some space, no need to keep the haplotagged bam file after phasing, if we are not using the bam later.
        keep = "--skip_final_phased_bam"
    shell:
        """
        run_pepper_margin_deepvariant call_variant \
            --bam {input.bam} \
            --fasta {input.genome} \
            --output_dir {params.output_dir} \
            --output_prefix {wildcards.sample} \
            --sample_name {wildcards.sample} \
            --phased_output {params.keep} \
            {params.mode} \
            --threads {threads} | tee -a {log}
        """