## use Delly to call somatic SV from ONT
## This workflow only take one pair of tumour-normal samples as input 
## you might consider building a multi-sample control panel with all normal samples if you have a cohort.

rule call_somatic_sv_delly:
    input:

    shell:
        """
        delly lr -y ont -o delly.bcf -g {input.genome} {input.tumour bam}
        """