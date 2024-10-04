import pandas as pd
from pathlib import Path
import sys

##TODO: pick the one from germline snv callers and somatic snv callers (some takes a long time, no need to run every tool), considering merging SVs results.

samples_df = pd.read_csv(config['samples'])
samples_df['sample_id'] = samples_df['sample_id'].astype(str) #sample id could be numbers only.
germline_samples_df = samples_df[samples_df['type']=='normal']

wildcard_constraints:
    sample="|".join(samples_df["sample_id"].unique()),
    run="|".join(samples_df["flowcell_id"].unique()),


## get all bam files of a sample, could be data from multiple flowcells/runs.
def get_bam_of_runs(wildcards):
    sample = wildcards.sample
    align_results_path = Path(f'analysis/bam/{sample}')
    flowcell_ids = set(samples_df[ (samples_df['sample_id'] == sample) ].flowcell_id)
    return [align_results_path / f'{run}.bam' for run in flowcell_ids]


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

    def get_clair3_output(self):
        return collect("analysis/snvs/clair3/{sample}/merge_output.vcf.gz", sample=self.germline_df['sample_id'].unique())

    def get_deepvariant_output(self):
        return collect("analysis/snvs/deepvariant/{sample}/{sample}.vcf.gz", sample=self.germline_df['sample_id'].unique())

    def get_deepsomatic_output(self):
        return [f"analysis/snvs/deepsomatic/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in self.pairs]

    def get_clairs_output(self):
        return [f"analysis/snvs/clairS/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in self.pairs]

    def get_snv_indel_calling_output(self):
        results = []
        for tool in self.config['snv_calling']['germline'] + self.config['snv_calling']['somatic']:
            if tool.lower() == 'deepvariant':
                results += self.get_deepvariant_output()
            elif tool.lower() == 'clair3':
                results += self.get_clair3_output()
            elif tool.lower() == 'deepsomatic':
                results += self.get_deepsomatic_output()
            elif tool.lower() == 'clairs':
                results += self.get_clairs_output()
            else:
                sys.exit(f"Error: snv_caller {tool} is not supported.")
        return results
    
    def get_nanomonsv_output(self):
        return [f"analysis/svs/nanomonsv/{pair['tumour']}.{pair['normal']}/{pair['tumour']}.{pair['normal']}.nanomonsv.sbnd.annot.proc.result.pass.txt" for pair in self.pairs]
    


def get_final_output():
    output_collector = OutputCollector(samples_df, config)
    aligned_bam = output_collector.get_alignment_output()
    snv = output_collector.get_snv_indel_calling_output()
    sv = output_collector.get_nanomonsv_output()
    return aligned_bam + snv + sv