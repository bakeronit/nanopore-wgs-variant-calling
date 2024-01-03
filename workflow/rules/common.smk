import pandas as pd
from pathlib import Path
import sys

samples_df = pd.read_csv(config['samples'])

wildcard_constraints:
    sample="|".join(samples_df["sample_id"])

MODE = config['basecalling_mode']

def get_basecalling_output():
    return expand(("analysis/ubam/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '/' + samples_df['flowcell_id'] + '.ubam'), mode=MODE)

def get_alignment_output():
    return expand(("analysis/bam/" + samples_df['flowcell_version'] + "/{mode}/" +  samples_df['sample_id'] + '.bam').unique(), mode=MODE)

def generate_paired_samples(df):
    """generate tumour-normal pairs for each donor"""
    donor_flowcell_df = df[['donor_id','flowcell_version']].drop_duplicates() # for each donor, and potentially each flowcell version
    pairs = []
    for index, row in donor_flowcell_df.iterrows():
        donor_id = row['donor_id']
        flowcell = row['flowcell_version']
        tumour_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'tumour')]['sample_id'].unique().tolist()
        normal_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'normal')]['sample_id'].unique().tolist()
        if len(tumour_sample) > 1 and len(normal_sample) > 1:
            raise ValueError('Multiple tumour vs multiple normal samples for one donor, could not determine the paired samples')
            sys.exit()
        pairs += [{'donor': donor_id, 'tumour': t, 'normal': n, 'flowcell_version': flowcell} for t in tumour_sample for n in normal_sample]
    return pairs

def get_snv_calling_output(tool):
    if tool in ['clairS','deepsomatic']:
        pairs = generate_paired_samples(samples_df)
        output = [f"analysis/snvs/{tool}/{pair['flowcell_version']}/{m}/{pair['tumour']}.{pair['normal']}/output.vcf.gz" for pair in pairs for m in MODE]
    if tool == 'clair3':
        clair3_path = Path('analysis/snvs/clair3')
        output = expand(( clair3_path/samples_df['flowcell_version']/'{mode}'/samples_df['sample_id']/'merge_output.vcf.gz').unique(), mode=MODE)
    if tool == 'deepvariant':
        deepvariant_path= Path('analysis/snvs/deepvariant/R10')  # deepvariant only works with R10 data.
        samples_r10_df = samples_df[samples_df['flowcell_version']=='R10']
        output = expand(( deepvariant_path/'{mode}'/samples_r10_df['sample_id']/(samples_r10_df['sample_id'] + ".vcf.gz")).unique(), mode=MODE)
        output += expand(( deepvariant_path/'{mode}'/samples_r10_df['sample_id']/(samples_r10_df['sample_id'] + ".g.vcf.gz")).unique(), mode=MODE)
    if tool == 'pepper':
        pepper_path = Path('analysis/snvs/pepper')
        output = expand(( pepper_path/samples_df['flowcell_version']/'{mode}'/samples_df['sample_id']/(samples_df['sample_id'] + '.phased.vcf.gz')).unique(), mode=MODE) 
    return output

def get_final_output():
    ## step1: gather data
    #final_output = ("raw_pod5/" + samples_df['flowcell'] + '/' + samples_df['sample_id'] + '/' + samples_df['run_id'] + '.done')
    
    ## step2: basecalling
    #basecalling_output = get_basecalling_output()

    ## step3: alignment
    #alignment_output = get_alignment_output()  

    ## step 3a: qc
    #final_output = expand(("analysis/qc/basecalling/" + samples_df['flowcell'] + "/{mode}/" +  samples_df['sample_id'] + '/' + samples_df['run_id'] + '.summary.txt'), mode=['hac','sup'])
    #final_output += expand(("analysis/qc/bam/" + samples_df['flowcell'] + "/{mode}/" +  samples_df['sample_id'] + '{suffix}'), mode=['hac','sup'], suffix=['.mosdepth.summary.txt','.mosdepth.global.dist.txt','.bamcov.txt']) 

    ## step4: snv calling
    snv_output = []
    for tool in ['clairS','deepsomatic','pepper','deepvariant']:
        snv_output += get_snv_calling_output(tool)
    
    ## methylation
    mod_output = expand( ('analysis/mod/' + samples_df['flowcell_version'] + '/{mode}/' + samples_df['sample_id'] + '.bed.gz').unique(), mode=MODE)
    mod_output += expand( ('analysis/mod/' + samples_df['flowcell_version'] + '/{mode}/' + samples_df['sample_id'] + '.mod_summary.txt').unique(), mode=MODE)

    ## benchmarks
    #final_output += expand(generate_clairS_paired_samples(samples_df, mode='{mode}', type = 'benchmark'), mode=['hac','sup'])
    #germline_df = samples_df[samples_df['type'] == 'normal']
    #final_output +=  expand(("analysis/benchmarks/snvs/germline/{caller}/" + germline_df['flowcell'] + "/{mode}/" + germline_df['sample_id'] + "/summary.txt").unique(),caller = ['clair3','pepper'], mode=['hac','sup'] )
    #final_output +=  expand("analysis/benchmarks/snvs/germline/deepvariant/R10/{mode}/COLO829_BL/summary.txt",mode=['hac','sup'])
    final_output = mod_output + snv_output   # + sv_output + phased_output + qc_output + benchmark_output

    step = config['step']
    if step == 'basecalling':
        return get_basecalling_output()
    elif step == 'alignment':
        return get_alignment_output()
    elif step == 'snv':
        return snv_output
    elif step == 'all':
        return final_output
    else:
        raise ValueError("step can only be basecalling, alignment, snv, and all.")
        return None

## get both cannonical and mod model files for basecalling
def get_basecalling_models(wildcards):
    if wildcards.flowcell == "R9":
        if wildcards.mode == "hac":
            cn_model = config['ont_model']['r9']['hac']['canonical']
            remora_model = config['ont_model']['r9']['hac']['remora']
        elif wildcards.mode == "sup":
            cn_model = config['ont_model']['r9']['sup']['canonical']
            remora_model = config['ont_model']['r9']['sup']['remora']
    elif wildcards.flowcell == "R10":
        if wildcards.mode == "hac":
            cn_model = config['ont_model']['r10']['hac']['canonical']
            remora_model = config['ont_model']['r10']['hac']['remora']
        elif wildcards.mode == "sup":
            cn_model = config['ont_model']['r10']['sup']['canonical']
            remora_model = config['ont_model']['r10']['sup']['remora']
    else:
        print("Error: flowcell version has to be R9 or R10. mode has to be hac or sup")
        exit(0)
    
    return {"cn_model": cn_model, "remora_model": remora_model}

## get the path to raw output from sequencer, either fast5 or pod5
def get_raw_path(wildcards):
    raw_path = list(samples_df[ (samples_df['flowcell_id'] == wildcards.run) ].raw_path)
    return raw_path
    
## get all bam files of a sample, could be data from multiple flowcells/runs.
def get_bam_of_runs(wildcards):
    sample = wildcards.sample
    flowcell = wildcards.flowcell
    mode = wildcards.mode

    align_results_path = Path(f'analysis/bam/{flowcell}/{mode}/{sample}')
    flowcell_ids = list(samples_df[ (samples_df['flowcell_version'] == flowcell) & (samples_df['sample_id'] == sample) ].flowcell_id)

    return [align_results_path / f'{id}.bam' for id in flowcell_ids]

## get model file for clair3 snp calling
def get_clair3_model(wildcards):
    if wildcards.flowcell == 'R9':
        clair3_model = config['clair3']['model_r9']
    elif wildcards.flowcell == 'R10':
        clair3_model = config['clair3']['model_r10']
    else:
        print("Error: flowcell version has to be R9 or R10.")
        exit(0)
    
    return clair3_model
    
def generate_clairS_paired_samples(df, mode, type = 'snv'):
    donor_flowcell_df = df[['donor_id','flowcell_version']].drop_duplicates()
    clairs_outputs = []
    for index, row in donor_flowcell_df.iterrows():
        donor_id = row['donor_id']
        flowcell = row['flowcell_version']

        tumour_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'tumour')]['sample_id'].unique().tolist()
        normal_sample = df[(df['donor_id'] == donor_id) & (df['flowcell_version'] == flowcell) & (df['type'] == 'normal')]['sample_id'].unique().tolist()

        if type == 'snv':
            clairs_outputs += [f"analysis/snvs/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/output.vcf.gz"]
        else: ## type = 'benchmark'
            clairs_outputs += [f"analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/fp_fn.vcf"]
            clairs_outputs += [f"analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/fp.vcf"]
            clairs_outputs += [f"analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/summary.txt"]
            clairs_outputs += [f"analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/tp.vcf"]
            clairs_outputs += [f"analysis/benchmarks/snvs/somatic/clairS/{flowcell}/{mode}/{tumour_sample[0]}.{normal_sample[0]}/fn.vcf"]

    return clairs_outputs