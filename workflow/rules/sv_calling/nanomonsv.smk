### use nanomonsv to call somatic SVs. 
## set --single_bnd --use_racon, so it can identify single breakend sv as well.
## use --qv20 preset parameters for nanopore Q20 chemistry. 
## set --min_indel_size 10 (default 50) to see if it can get smallers sv.maybe not worth it.
## nanomonsv_filter_simple_repeat_svtype: Post filtering of simple repeat and add svtype to the txt results
## nanomonsv_postprocess_sbnd: annotate svs and classify sbnd.
###########

misc = config['nanomonsv']['misc_scripts_path']
rule nanomonsv_postprocess_sbnd:
    input:
        genome = config['reference']['file'],
        result = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.txt",
        sbnd_result = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd.result.txt",
        simple_repeat = config['nanomonsv']['simple_repeat']
    output:
        directory("analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd_vis"),
        annot = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.annot.proc.result.txt",
        sbnd_annot = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd.annot.proc.result.txt",
        annot_pass = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.annot.proc.result.pass.txt",
        sbnd_annot_pass = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd.annot.proc.result.pass.txt"
    params:
        prefix = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}"
    envmodules:
        "nanomonsv/0.7.1",    
        "R/4.3.1"  # need tidyverse and ggrepel
    shell:
        """
        python3 {misc}/subscript_postprocess_sbnd/integrate_sbnd.py {params.prefix} {input.genome} 
        python3 {misc}/add_simple_repeat.py {params.prefix}.nanomonsv.proc.result.txt {params.prefix}.nanomonsv.annot.proc.result.txt.tmp {input.simple_repeat}
        nanomonsv insert_classify {params.prefix}.nanomonsv.annot.proc.result.txt.tmp {output.annot} {input.genome} --genome_id hg38
        python3 {misc}/subscript_postprocess_sbnd/add_simple_repeat_sbnd.py {params.prefix}.nanomonsv.sbnd.proc.result.txt {input.simple_repeat} > {output.sbnd_annot}
        Rscript {misc}/subscript_postprocess_sbnd/plot_sbnd_vis.R {params.prefix} {params.prefix}.nanomonsv.sbnd_vis
        rm -rf {params.prefix}.nanomonsv.annot.proc.result.txt.tmp

        head -n 1 {output.annot} > {output.annot_pass}
        tail -n +2 {output.annot} | grep PASS >> {output.annot_pass}
        head -n 1 {output.sbnd_annot} > {output.sbnd_annot_pass}
        tail -n +2 {output.sbnd_annot} | grep PASS >> {output.sbnd_annot_pass}
        """

rule nanomonsv_filter_simple_repeat_svtype:
    input:
        result = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.txt",
        simple_repeat = config['nanomonsv']['simple_repeat']
    output:
        filt = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.txt",
        passed = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.pass.txt",
        svtype = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.filt.pass.svtype.txt"
    threads: 1
    resources:
        mem = 2,
        walltime = 1
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        python3 {misc}/add_simple_repeat.py {input.result} {output.filt} {input.simple_repeat}

        head -n 1 {output.filt} > {output.passed} 
        tail -n +2 {output.filt} |grep PASS >> {output.passed}

        python3 {misc}/sv_type.py {output.passed} {output.svtype}
        """

rule call_somatic_sv_nanomonsv_get:
    input:
        [f"analysis/svs/nanomonsv/{{flowcell}}/{{mode}}/{{sample_t}}/{{sample_t}}.{type}.sorted.bed.gz" for type in ['bp_info','deletion','insertion']],
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}/{sample_t}.rearrangement.sorted.bedbp.gz",
        [f"analysis/svs/nanomonsv/{{flowcell}}/{{mode}}/{{sample_n}}/{{sample_n}}.{type}.sorted.bed.gz" for type in ['bp_info','deletion','insertion']],
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_n}/{sample_n}.rearrangement.sorted.bedbp.gz",
        genome = config['reference']['file'],
        tumour_bam = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam",
        tumour_bai = "analysis/bam/{flowcell}/{mode}/{sample_t}.bam.bai",
        normal_bam = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        normal_bai = "analysis/bam/{flowcell}/{mode}/{sample_n}.bam",
        control_panel_path = config['nanomonsv']['control_panel_path'],
    output:
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.result.txt",
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}.{sample_n}/{sample_t}.{sample_n}.nanomonsv.sbnd.result.txt"
    params:
        tumour_prefix = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_t}/{sample_t}",
        normal_prefix = "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample_n}/{sample_n}",
        panel_prefix="hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control"
    threads: 20
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        nanomonsv get {params.tumour_prefix} {input.tumour_bam} {input.genome} \
        --control_prefix {params.normal_prefix} --control_bam {input.normal_bam} \
        --single_bnd --use_racon --min_indel_size 10 --qv20 \
        --control_panel_prefix {input.control_panel_path}/{params.panel_prefix} --processes {threads} 
        """

rule nanomonsv_parse:
    input:
        bam = "analysis/bam/{flowcell}/{mode}/{sample}.bam",
        bai="analysis/bam/{flowcell}/{mode}/{sample}.bam.bai",
        #genome = config['reference']['file']   # nnmsv supports cram from v0.7.0. for cram input, --reference_fasta is recommended.
    output:
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}.bp_info.sorted.bed.gz",
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}.bp_info.sorted.bed.gz.tbi",
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}.deletion.sorted.bed.gz",
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}.insertion.sorted.bed.gz",
        "analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}.rearrangement.sorted.bedpe.gz",
    params:
        prefix="analysis/svs/nanomonsv/{flowcell}/{mode}/{sample}/{sample}"
    threads: 1
    resources:
        mem = 24,
        walltime = 48
    envmodules:
        "nanomonsv/0.7.1"
    shell:
        """
        nanomonsv parse {input.bam} {params.prefix}
        """