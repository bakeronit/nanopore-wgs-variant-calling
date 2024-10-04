# trim the heading of 40bp and tailing of 20bp for each read.
# after trimmming, need to fix the MM tag.
modkit = config['modkit']

rule trim_reads:
    input:
        "analysis/ubam/{sample}/{run}.ubam"
    output:
        temp("analysis/ubam/{sample}/{run}.trimmed.sorted.ubam")
    threads: 12
    resources:
        mem = 20,
        walltime = 24
    envmodules:
        "samtools/1.17"
    params:
        headcrop = config['filter']['headcrop'],
        tailcrop = config['filter']['tailcrop'],
        minlen = config['filter']['minlen'],
        trim_len = config['filter']['headcrop'] + config['filter']['tailcrop']
    shell:
        """
        samtools fastq -T"*" {input} | \
        chopper --headcrop {params.headcrop} --tailcrop {params.tailcrop} --minlength {params.minlen} --threads {threads} | \
        awk '/^@.*MN:i:[0-9]+/ {{ match($0, /MN:i:([0-9]+)/, num);sub(/MN:i:[0-9]+/, "MN:i:" num[1]-{params.trim_len}) }}  1' |\
        samtools import -T"*" - | \
        samtools sort -@{threads} -n > {output}
        """

rule sorted_original_ubam:
    input:
        "analysis/ubam/{sample}/{run}.ubam"
    output:
        temp("analysis/ubam/{sample}/{run}.sorted.ubam")
    threads: 12
    resources:
        mem = 20,
        walltime = 12
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools sort -@{threads} -n {input} > {output}
        """

rule repair_MMtag:
    input:
        original = "analysis/ubam/{sample}/{run}.sorted.ubam",
        trimmed = "analysis/ubam/{sample}/{run}.trimmed.sorted.ubam"
    output:
        temp("analysis/ubam/{sample}/{run}.trimmed_repaired.ubam")
    log:
        "logs/modkit_repair/{sample}.{run}.log"
    resources:
        mem = 20,
        walltime = 24
    threads: 10
    retries: 3
    shell:
        """
        {modkit} repair --donor-bam {input.original} --acceptor-bam {input.trimmed} \
            --output-bam {output} --log-filepath {log} --threads {threads}
        """
