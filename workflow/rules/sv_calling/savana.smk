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
    threads: 30
    resources:
        mem = 600,
        walltime = 24
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda  activate ~/working/local/micromanba_envs/savana1.20
        set -eu
        
        if [ "$(ls -A {params.outdir})" ]; then   # need to make sure the outdir is empty
            rm -f {params.outdir}/*
        fi

        savana run --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --length {params.min_svlen} \
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
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda  activate ~/working/local/micromanba_envs/savana1.20
        set -eu

        savana classify \
        --vcf {input.vcf} \
        --output {output.vcf} \
        --somatic_output {output.somatic_vcf} \
        --ont \
        --cna_rescue \
        --threads {threads} &> {log}
        """

rule call_somatic_sv_savana_cna:
    input:
        phased_vcf = "analysis/snvs/clair3/{sample_n}/phased_merge_output.vcf.gz",
        tumour_bam="analysis/bam/{sample_t}.bam",
        tumour_bai="analysis/bam/{sample_t}.bam.bai",
        normal_bam="analysis/bam/{sample_n}.bam",
        normal_bai="analysis/bam/{sample_n}.bam.bai",
        genome=config['reference']['file'],
        blacklist=config['blacklist'],
    output:
        "analysis/cnvs/savana/{sample_t}.{sample_n}/10kbp_bin_ref_all_{sample_t}.{sample_n}withVCF.bed",
        multiext("analysis/cnvs/savana/{sample_t}.{sample_n}/{sample_t}.{sample_n}", \
        "_allele_counts_hetSNPs.bed", \
        "_fitted_purity_ploidy.tsv", \
        "_phased_het_snps.bed", \
        "_ranked_solutions.tsv", \
        "_read_counts_filtered.tsv", \
        "_read_counts_mnorm_log2r_segmented.tsv", \
        "_read_counts_mnorm_log2r_smoothened_sl10_t0.025.tsv", \
        "_read_counts_mnorm_log2r.tsv", \
        "_read_counts.tsv", \
        "_segmented_absolute_copy_number.tsv" )
    params:
        outdir="analysis/cnvs/savana/{sample_t}.{sample_n}",
    log:
        "logs/savana/{sample_t}.{sample_n}.cna.log"
    benchmark:
        "benchmarks/savana/{sample_t}.{sample_n}.cna.benchmark.txt"
    threads: 30
    resources:
        mem = 400,
        walltime = 48
    envmodules:
        "conda-envs/base"
    shell:
        """
        set +eu
        conda activate ~/working/local/micromanba_envs/savana1.20
        set -eu

        savana cna --tumour {input.tumour_bam} \
        --normal {input.normal_bam} \
        --ref {input.genome} \
        --sample {wildcards.sample_t}.{wildcards.sample_n} \
        --phased_vcf {input.phased_vcf} \
        --breakpoints {input.phased_vcf} \
        --blacklist {input.blacklist} \
        --threads {threads} \
        --outdir {params.outdir} &> {log}
        """
        
    