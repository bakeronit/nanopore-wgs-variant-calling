rule call_germline_snv_deepvariant_vcf_stats_report:
    input:
        "analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz",
    output:
        "analysis/snvs/deepvariant/{sample}/{sample}.visual_report.html"
    params:
        outbase = "analysis/snvs/deepvariant/{sample}/{sample}"
    threads: 1
    resources:
        mem = 10,
        walltime = 10
    container: DEEPVARIANT_CPU_sif
    shell:
        """
        vcf_stats_report \
        --input_vcf {input} \
        --outfile_base {params.outbase}
        """

rule filter_passed_variants:
    input:
        "analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz"
    output:
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.passed.vcf.gz",
        tbi = "analysis/snvs/deepvariant/{sample}/{sample}.passed.vcf.gz.tbi"
    threads: 1
    resources:
        mem = 6,
        walltime = 2
    shell:
        """
        bcftools view -i 'FILTER="PASS"' -Oz -o {output.vcf} {input}
        bcftools index -t {output.vcf}
        """
    
#rule get_chromosomes:
#    input:
#        "analysis/snvs/deepvariant/{sample}/{sample}.passed.vcf.gz"
#    output:
#        pipe("analysis/snvs/deepvariant/{sample}/chromosomes.txt")
#    threads: 1
#    resources:
#        mem = 6,
#        walltime = 2
#    shell:
#        """
#        bcftools index --stats {input} | cut -f 1 > {output}
#        """

rule whatshap_phasing:
    input:
        reference = config['reference']['file'],
        bam = "analysis/bam/{sample}.bam",
        vcf = "analysis/snvs/deepvariant/{sample}/{sample}.passed.vcf.gz"
    output:
        "analysis/snvs/deepvariant/{sample}/{sample}.passed.phased.vcf.gz"
    threads: 12
    resources:
        mem = 24,
        walltime = 24
    shell:
        """
        whatshap phase \
        --output {output} \
        --reference {input.reference} \
        --ignore-read-groups \
        {input.vcf} \
        {input.bam}
        """

rule index_phased_vcf:
    input:
       rules.whatshap_phasing.output
    output:
        "analysis/snvs/deepvariant/{sample}/{sample}.passed.phased.vcf.gz.tbi"
    threads: 1
    resources:
        mem = 1,
        walltime = 1
    shell:
        """
        bcftools index -t {input}
        """