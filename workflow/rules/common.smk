import sys
import pandas as pd
from pathlib import Path

container: "/mnt/backedup/home/jiaZ/working/containers/definitions/long_read_wgs_pipeline.sif"

samples_df = pd.read_csv(config['samples'])
samples_df['sample_id'] = samples_df['sample_id'].astype(str) #sample id could be numbers only.

wildcard_constraints:
    sample="|".join(samples_df["sample_id"].unique()),
    run="|".join(samples_df["flowcell_id"].unique()),

if config['run_mode'] in ['somatic', 'all']:
    wildcard_constraints:
        sample_t="|".join(samples_df[samples_df['type']=='tumour']['sample_id'].unique()),
        sample_n="|".join(samples_df[samples_df['type']=='normal']['sample_id'].unique())

def get_bam_of_flowcells(wildcards):
    sample = wildcards.sample
    align_results_path = Path(f'analysis/bam/{sample}')
    flowcell_ids = set(samples_df[ (samples_df['sample_id'] == sample) ].flowcell_id)
    return [align_results_path / f'{flowcell}.bam' for flowcell in flowcell_ids]

def get_phased_vcf(wildcards):
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards.sample_t
    donor_id = samples_df[samples_df['sample_id']==sample]['donor_id'].iloc[0]
    normal_sample_id = samples_df[(samples_df['donor_id']==donor_id) & (samples_df['type']=='normal')]['sample_id'].iloc[0]
    phased_snv_from = config['snv_calling']['phased_snv_from']
    if phased_snv_from not in config['snv_calling']['germline']:
        print("Warning: phased vcf is not generated from the germline snv caller you used for this workflow.")
        phased_snv_from = config['snv_calling']['germline'][0]
    match phased_snv_from.lower():
        case 'clair3':
            return f"analysis/snvs/clair3/{normal_sample_id}/phased_merge_output.vcf.gz"
        case 'pepper':
            return f"analysis/snvs/pepper/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case 'deepvariant':
            return f"analysis/snvs/deepvariant/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case _:
            raise ValueError(f"{phased_snv_from} is not supported!")

def get_alignment_output(df):
    return collect("analysis/bam/{sample}.bam", sample=df.sample_id.unique())

def get_qc_output(df):
    return collect("analysis/qc/bam/{sample}.mosdepth.global.dist.txt", sample=df['sample_id'].unique()) + \
        collect("analysis/qc/bam/{sample}.mosdepth.summary.txt", sample=df['sample_id'].unique()) + \
        collect("analysis/qc/bam/{sample}.bamcov.txt", sample=df['sample_id'].unique())

def generate_paired_samples(df):
    """
    Generate test-control pairs for each donor, there might be multiple test samples with one control sample.
    For cases with multiple control samples, it will raise an error.
    """
    pairs = []
    if config['run_mode'] not in ['somatic', 'all']: 
        return pairs

    for donor_id in df['donor_id'].unique():
        tumour_sample = df[(df['donor_id'] == donor_id) & (df['type'] == 'tumour')]['sample_id'].unique()
        normal_sample = df[(df['donor_id'] == donor_id) & (df['type'] == 'normal')]['sample_id'].unique()
        if len(normal_sample) > 1:
            raise ValueError(f'Multiple control samples for {donor_id}, could not determine the paired samples, please check your sample sheet.')
        pairs += [{'donor': donor_id, 'tumour': t, 'normal': n} for t in tumour_sample for n in normal_sample]
    return pairs

def get_snv_indel_output(df, caller):
    pairs = generate_paired_samples(df)
    results = {
        'pepper': collect("analysis/snvs/pepper/{sample}/{sample}.vcf.gz", sample=df.sample_id.unique()),
        'clair3': collect("analysis/snvs/clair3/{sample}/merge_output.vcf.gz", sample=df.sample_id.unique()),
        'deepvariant': collect("analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz", sample=df.sample_id.unique()),
        'clairs': [f"analysis/snvs/clairS/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in pairs],
        'deepsomatic': [f"analysis/snvs/deepsomatic/{pair['tumour']}.{pair['normal']}/output.somatic.vcf.gz" for pair in pairs]
    }
    if results.get(caller.lower()) is None:
        raise ValueError(f"Error: SNV caller {caller} is not supported.")
    return results[caller.lower()]

def get_sv_output(df, caller):
    pairs = generate_paired_samples(df)
    results = {
        'sniffles': collect("analysis/svs/sniffles/{sample}/{sample}.vcf", sample=df.sample_id.unique()), 
        'nanomonsv': [f"analysis/svs/nanomonsv/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.nanomonsv.annot.proc.result.txt" for pair in pairs],
        'severus': [f"analysis/svs/severus/{pair['tumour']}.{pair['normal']}/somatic_SVs/severus_somatic.vcf" for pair in pairs],
        'savana': [f"analysis/svs/savana/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.classified.somatic.vcf" for pair in pairs]
    }
    if results.get(caller.lower()) is None:
        raise ValueError(f"Error: SV caller {caller} is not supported.")
    return results[caller.lower()]

def get_final_output():
    run_mode = config['run_mode']
    final_results = [
        get_alignment_output(samples_df),
        get_qc_output(samples_df)
    ]
    if run_mode in ['germline', 'all']:
        final_results += [get_snv_indel_output(samples_df, caller) for caller in config['snv_calling']['germline']]
        final_results += [get_sv_output(samples_df, caller) for caller in config['sv_calling']['germline']]
    if run_mode in ['somatic', 'all']:
        final_results += [get_sv_output(samples_df, caller) for caller in config['sv_calling']['somatic']]
    return final_results