## call somatic SV with long-read data using SAVANA. https://github.com/cortes-ciriano-lab/savana
## set --length: Minimum length SV to consider (default=30),
## set --mapq: Minimum MAPQ of reads to consider (default=5),
## set --depth: Minumum number of supporting reads from tumour OR normal to consider variant (default=3)
## set --buffer: Buffer to add when clustering adjacent (non-insertion) potential breakpoints (default=10)
#####

rule call_somatic_sv_savana_run:
    input:
        tumour_bam="analysis/bam/{sample_t}.bam",
        tumour_bai="analysis/bam/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{sample_n}.bam",
        normal_bai="analysis/bam/{sample_n}.bam.bai",
        genome=config['reference']['file'],
    output:
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.bedpe",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints_read_support.tsv",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.inserted_sequences.fa",
    params:
        outdir="analysis/svs/savana/{sample_t}.{sample_n}",
        min_svlen=50
    log:
        "logs/savana/{sample_t}.{sample_n}.run.log"
    benchmark:
        "benchmarks/savana/{sample_t}.{sample_n}.run.benchmark.txt"
    threads: 24 # to run three jobs a time in bigmem node
    resources:
        mem = 600,
        walltime = 24
    shell:
        """
        if [ "$(ls -A {params.outdir})" ]; then   # need to make sure the outdir is empty
            rm -f {params.outdir}/*
        fi

        savana run --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --length {params.min_svlen} \
        --single_bnd \
        --threads {threads} \
        --outdir {params.outdir} &> {log}
        """

rule call_somatic_sv_savana_classify:
    input:
        vcf = "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
    output:
        vcf = "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.sv_breakpoints.vcf",
        somatic_vcf = "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.classified.somatic.vcf"
    log:
        "logs/savana/{sample_t}.{sample_n}.classify.log"
    benchmark:
        "benchmarks/savana/{sample_t}.{sample_n}.classify.benchmark.txt"
    params:
        outdir="analysis/svs/savana/{sample_t}.{sample_n}"
    threads: 12
    resources:
        mem = 20,
        walltime = 48
    shell:
        """
        savana classify \
        --vcf {input.vcf} \
        --output {output.vcf} \
        --somatic_output {output.somatic_vcf} \
        --ont \
        --cna_rescue \
        --threads {threads} &> {log}
        """

# call somatic CNV with savana cna subcommand
# the dependency of phased vcf file is not needed since version 1.3.0
# the output files names are also changed, need to keep an eye on updates.
rule call_somatic_cnv_savana_cna:
    input:
        phased_vcf = get_phased_vcf,
        breakpoints = "analysis/svs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}.sv_breakpoints.vcf",
        tumour_bam="analysis/bam/{sample_t}.bam",
        tumour_bai="analysis/bam/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{sample_n}.bam",
        normal_bai="analysis/bam/{sample_n}.bam.bai",
        genome=config['reference']['file'],
        blacklist=config['annotation']['blacklist'],
    output:
        "analysis/cnvs/savana/{sample_t}.{sample_n}/10kbp_bin_ref_all_{sample_t}.{sample_n}_with_SV_breakpoints.bed",
        multiext("analysis/cnvs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}", \
        "_allele_counts_hetSNPs.bed", \
        "_fitted_purity_ploidy.tsv", \
        "_ranked_solutions.tsv", \
        "_read_counts_mnorm_log2r_segmented.tsv", \
        "_segmented_absolute_copy_number.tsv" )
    params:
        outdir="analysis/cnvs/savana/{sample_t}.{sample_n}",
    log:
        "logs/savana/{sample_t}.{sample_n}.cna.log"
    benchmark:
        "benchmarks/savana/{sample_t}.{sample_n}.cna.benchmark.txt"
    threads: 24
    resources:
        mem = 200,
        walltime = 48
    shell:
        """
        if [ "$(ls -A {params.outdir})" ]; then   
            rm -f {params.outdir}/*
        fi

        savana cna \
        --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --snp_vcf {input.phased_vcf} \
        --breakpoints {input.breakpoints} \
        --blacklist {input.blacklist} \
        --threads {threads} \
        --outdir {params.outdir} &> {log}
        """
        
    