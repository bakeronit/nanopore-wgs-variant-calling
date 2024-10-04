import pandas as pd
from pathlib import Path
import sys

##TODO: pick the one from germline snv callers and somatic snv callers (some takes a long time, no need to run every tool), considering merging SVs results.

samples_df = pd.read_csv(config['samples'])
samples_df['sample_id'] = samples_df['sample_id'].astype(str) #sample id could be numbers only.
wildcard_constraints:
    sample="|".join(samples_df["sample_id"].unique()),
    run="|".join(samples_df["flowcell_id"].unique()),

def get_alignment_output():
    return collect("analysis/bam/{sample}.bam", sample=samples_df['sample_id'].unique())

## get all bam files of a sample, could be data from multiple flowcells/runs.
def get_bam_of_runs(wildcards):
    sample = wildcards.sample
    align_results_path = Path(f'analysis/bam/{sample}')
    flowcell_ids = set(samples_df[ (samples_df['sample_id'] == sample) ].flowcell_id)
    return [align_results_path / f'{flowcell}.bam' for flowcell in flowcell_ids]

def get_clair3_output():
    return collect("analysis/snvs/clair3/{sample}/merge_output.vcf.gz", sample=samples_df['sample_id'].unique())

def get_deepvariant_output():
    return collect("analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz", sample=samples_df['sample_id'].unique())

def get_germline_snv_indel_calling_output():
    results = []
    for tool in config['snv_calling']['germline']:
        if tool == 'deepvariant':
            results += get_deepvariant_output()
        elif tool == 'clair3':
            results += get_clair3_output()
        else:
            sys.exit(f"Error: snv_caller {tool} is not supported.")
    return results

def get_final_output():
    aligned_bam = get_alignment_output()
    germline_snv = get_germline_snv_indel_calling_output()
    
    return aligned_bam + germline_snv