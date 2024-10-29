import pandas as pd
from pathlib import Path
import sys

samples_df = pd.read_csv(config['samples'])
samples_df['sample_id'] = samples_df['sample_id'].astype(str) #sample id could be numbers only.

wildcard_constraints:
    sample="|".join(samples_df["sample_id"].unique()),
    run="|".join(samples_df["flowcell_id"].unique()),
    sample_t="|".join(samples_df[samples_df['type']=='tumour']['sample_id'].unique()),
    sample_n="|".join(samples_df[samples_df['type']=='normal']['sample_id'].unique())

## get all bam files of a sample
def get_bam_of_runs(wildcards):
    sample = wildcards.sample
    align_results_path = Path(f'analysis/bam/{sample}')
    flowcell_ids = set(samples_df[ (samples_df['sample_id'] == sample) ].flowcell_id)
    return [align_results_path / f'{run}.bam' for run in flowcell_ids]

def get_phased_vcf(wildcards):
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards.sample_t
    donor_id = samples_df[samples_df['sample_id']==sample].donor_id.tolist()[0]
    normal_sample_id = samples_df[(samples_df['donor_id']==donor_id) & (samples_df['type']=='normal')].sample_id.tolist()[0]
    if config['phased_snv_from'] not in config['snv_calling']['germline']:
        print("Warning: phased vcf is not generated from the germline snv caller you used for this workflow.")
    match config['phased_snv_from'].lower():
        case 'clair3':
            return f"analysis/snvs/clair3/{normal_sample_id}/phased_merge_output.vcf.gz"
        case 'pepper':
            return f"analysis/snvs/pepper/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case 'deepvariant':
            return f"analysis/snvs/deepvariant/{normal_sample_id}/{normal_sample_id}.phased.vcf.gz"
        case _:
            raise ValueError("haplotagged bam should only be generated from clair3, pepper, or deepvariant!")


class OutputCollector:
    def __init__(self, df: pd.DataFrame, config):
        self.df = df
        self.config = config
        self.germline_df = df[df['type'] == 'normal']
        self.pairs = self.generate_paired_samples()
    
    def get_alignment_output(self):
        return collect("analysis/bam/{sample}.bam", sample=self.df['sample_id'].unique())

    def generate_paired_samples(self):
        """generate tumour-normal pairs for each donor, there might be multiple tumour for one normal sample"""
        pairs = []
        for donor_id in self.df['donor_id'].unique():
            tumour_sample = self.df[(self.df['donor_id'] == donor_id) & (self.df['type'] == 'tumour')]['sample_id'].unique().tolist()
            normal_sample = self.df[(self.df['donor_id'] == donor_id) & (self.df['type'] == 'normal')]['sample_id'].unique().tolist()
            if len(tumour_sample) > 1 and len(normal_sample) > 1:
                raise ValueError('Multiple tumour vs multiple normal samples for one donor, could not determine the paired samples')
            pairs += [{'donor': donor_id, 'tumour': t, 'normal': n} for t in tumour_sample for n in normal_sample]
        return pairs

    def get_snv_indel_calling_output(self):
        results = []
        for tool in self.config['snv_calling']['germline'] + self.config['snv_calling']['somatic']:
            if tool.lower() == 'deepvariant':
                results += collect("analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz", sample=self.germline_df['sample_id'].unique())
            elif tool.lower() == 'clair3':
                results += collect("analysis/snvs/clair3/{sample}/merge_output.vcf.gz", sample=self.germline_df['sample_id'].unique())
            elif tool.lower() == 'deepsomatic':
                results += [f"analysis/snvs/deepsomatic/{pair['tumour']}.{pair['normal']}/output.somatic.vcf.gz" for pair in self.pairs]
            elif tool.lower() == 'clairs':
                results += [f"analysis/snvs/clairS/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in self.pairs]
            else:
                sys.exit(f"Error: snv_caller {tool} is not supported.")
        return results
    
    def get_savana_output(self, type='sv'):
        if type == 'sv':
            return [f"analysis/svs/savana/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.classified.somatic.vcf" for pair in self.pairs]
        elif type in ['cnv', 'cna']:
            return [f"analysis/cnvs/savana/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}_segmented_absolute_copy_number.tsv" for pair in self.pairs]
        else:
            raise ValueError(f"Type {type} is not supported.")

    def get_sv_calling_output(self):
        nanomonsv = [f"analysis/svs/nanomonsv/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.nanomonsv.sbnd.annot.proc.result.pass.txt" for pair in self.pairs]
        savana = self.get_savana_output(type='sv')
        severus = [f"analysis/svs/severus/{pair['tumour']}.{pair['normal']}/somatic_SVs/severus_somatic.vcf" for pair in self.pairs]

        return nanomonsv + savana + severus
    
    def get_qc_output(self):
        return collect("analysis/qc/bam/{sample}.mosdepth.global.dist.txt", sample=self.df['sample_id'].unique()) + \
            collect("analysis/qc/bam/{sample}.mosdepth.summary.txt", sample=self.df['sample_id'].unique()) + \
            collect("analysis/qc/bam/{sample}.bamcov.txt", sample=self.df['sample_id'].unique())  

def get_final_output(step='all'):
    output_collector = OutputCollector(samples_df, config)
    aligned_bam = output_collector.get_alignment_output()
    snv = output_collector.get_snv_indel_calling_output()
    sv = output_collector.get_sv_calling_output()
    cnv = output_collector.get_savana_output(type="cnv")
    qc = output_collector.get_qc_output()
    final_output = []
    match step:
        case "align":
            final_output = aligned_bam + qc
        case "snv":
            final_output = snv
        case "sv":
            final_output = sv
        case "all":
            final_output = snv + sv + qc #+ cnv JZ: can't test cnv due to small dataset.
        case _:
            raise ValueError(f"Step {step} is not supported.")

    return final_output