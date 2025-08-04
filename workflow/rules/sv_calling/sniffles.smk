rule call_germline_sv_sniffles:
    input:
        ref = config['reference']['file'],
        bam = "analysis/bam/{sample}.bam"
    output:
        "analysis/svs/sniffles/{sample}/{sample}.vcf"
    params:
        minsvlen = config['params']['minsvlen']
    log:
        "logs/sniffles/{sample}.log"
    benchmark:
        "benchmarks/sniffles/{sample}.benchmark.txt"
    threads: 24
    resources:
        mem = 48,
        walltime = 12
    shell:
        """
        sniffles --input {input.bam} \
            --vcf {output} \
            --reference {input.ref} \
            --minsvlen {params.minsvlen} \
            --allow-overwrite \
            --threads {threads} &> {log}
        """