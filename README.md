# Nanopore whole-genome sequencing variant calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.20.1-brightgreen.svg)](https://snakemake.github.io) [![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3100/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Sylabs SIF](https://img.shields.io/badge/Container-Sylabs%20SIF-5E72E4?style=flat-square)](https://cloud.sylabs.io/library/jiazhang/workflows/long_read_wgs_pipeline)

A Snakemake workflow for germline and somatic variant calling from Oxford Nanopore long-read whole-genome sequencing data of paired tumour/normal samples.

## Table of Contents

- [Overview](#overview)
- [Installation and Setup](#installation-and-setup)
  - [Requirements](#requirements)
  - [Container Setup](#container-setup)
  - [Download Required Models and Annotations](#download-required-models-and-annotations)
    - [Reference Genome](#reference-genome)
    - [Tool-Specific Models](#tool-specific-models)
    - [Tool-Specific Annotation Files](#tool-specific-annotation-files)
- [Configuration](#configuration)
  - [Sample Sheets](#sample-sheets)
  - [Workflow Configuration](#workflow-configuration)
- [Usage](#usage)
  - [Running the Complete Workflow](#running-the-complete-workflow)
  - [Monitoring Progress](#monitoring-progress-optional)
- [Output Structure](#output-structure)
  - [Key Output Files](#key-output-files)
- [Workflow Components](#workflow-components)
  - [Supported Tools](#supported-tools)
  - [Analysis Modes](#analysis-modes)
- [Resource Requirements](#resource-requirements)
  - [Compute Resources](#compute-resources)
  - [Optimising Workflow Resources](#optimising-the-workflow-resources-requsts-on-hpc)
- [Citation](#citation)
- [License](#license)

## Overview

![Workflow rulegraph](docs/rulegraph.png?width=600px)

This workflow processes Oxford Nanopore long-read WGS sequencing data to identify:
- **Small variants**: SNVs and indels in tumour and normal tissue
- **Structural variants**: Structural variants in tumour and normal tissue
- **Copy number alterations**: Segmental changes in copy number across the genome in tumour tissue
- **Phased variants**: Haplotype information for germline variants
- **Methylation**: DNA methylation at CpG sites in tumour and normal tissue

## Installation and Setup

### Requirements
- **Snakemake ≥8.20.1**
- **Singularity/Apptainer** (for containerized tools on HPC)
- **Python 3.11+**
- **Pandas**

### Container Setup

#### Option 1: Use Pre-built Container (Recommended)

```bash
apptainer pull library://jiazhang/workflows/long_read_wgs_pipeline:1.1
```

#### Option 2: Build Container Locally if new versions are needed

With [environment.yaml](workflow/envs/environment.yaml) and [long_read_wgs.def](workflow/envs/long_read_wgs.def) you can build the container locally.

```bash
# Build the container from definition file
apptainer build long_read_wgs_pipeline.sif long_read_wgs.def
```

### Download Required Models and Annotations

#### Reference Genome

A constant reference genome and indexed genome is required for all samples. Note that this pipeline is highly adapted to **GRCh38**.

```bash
# Download and index GRCh38 reference
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference.fasta

# Index reference
minimap2 -d reference.fasta.mmi reference.fasta
samtools faidx reference.fasta
```

#### Tool-Specific Models

- Nanopore rerio provides up-to-date pre-trained models for [Clair3](https://github.com/nanoporetech/rerio/tree/master/clair3_models).
- [ClairS](https://github.com/HKU-BAL/ClairS) `--platform` parameter changes with ONT models, need to keep it up-to-date. 

These four tools are not included in the pipeline container, but each tool has a containerized version available:

 - [Clair3](https://hub.docker.com/r/hkubal/clair3/tags)
 - [ClairS](https://hub.docker.com/r/hkubal/clairs/tags)
 - [DeepSomatic](https://hub.docker.com/r/google/deepvariant/tags)
 - [DeepVariant](https://hub.docker.com/r/google/deepvariant/tags)

#### Tool-Specific Annotation Files

**Nanomonsv** 
 - A panel SVs of control ([Pon](https://zenodo.org/record/7017953)) contains 30 Nanopore sequencing data from the Human Pangenome Reference Consortium. 
 - Pre-built [simple repeats](https://github.com/friend1ws/nanomonsv/tree/master/resource/simple_repeats) bed file.
 - Pre-built [LINE1_db](https://github.com/friend1ws/nanomonsv/tree/master/resource/LINE1_db)
 - Post-processing scripts: https://github.com/friend1ws/nanomonsv/tree/master/misc

**Severus VNTRs**

 - Variable number tandem repeats: [VNTRs](https://github.com/KolmogorovLab/Severus/tree/main/vntrs) annotations for common references, custom annotations can be created using [findTandemRepeats](https://github.com/PacificBiosciences/pbsv/blob/master/annotations/findTandemRepeats)
 - [PON](https://github.com/KolmogorovLab/Severus/tree/main/pon), used for **tumour only** cases. 

 **Blacklists**
 The ENCODE Blacklist, this is used in Delly and SAVANA for excluding problematic regions of the reference genome. 
 - [Blacklists](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)

## Configuration

### Sample Sheets

Create `config/samples.csv` with your sample information:
```csv
sample_id,flowcell_id,flowcell_version,donor_id,type
sample_t,PAQ23729,R10,PATIENT01,tumour
sample_n,PAQ24563,R10,PATIENT01,normal
```

> For QIMR Promethion sequencing, this sample file can be generated as in [nanopore basecalling pipeline](https://github.com/bakeronit/nanopore_dorado_basecalling), and manually add the type column indicating which one is the test (tumour) and which one is the control (normal).

**Required columns:**
- `sample_id`: Unique identifier for each sample
- `flowcell_id`: Flowcell barcode (e.g., PAQ23729)
- `flowcell_version`: Flowcell chemistry (R9/R10, best supporting R10)
- `donor_id`: Patient/donor identifier (must match between tumour/normal pairs)
- `type`: Sample type (`normal` or `tumour`; if running in germline mode, samples are processed as single samples with no somatic calling)

### Workflow Configuration

Edit `config/config.yaml`:
```yaml
# Sample and analysis configuration
samples: /path/to/samples.csv
run_mode: all
container: /path/to/long_read_wgs_pipeline.sif # path to container or library url
pull_containers: false # if running on qimrb cluster, no need to pull containers
basecalling_dir: /path/to/work_dir/to/basecalling_pipeline

snv_calling:
    phased_snv_from: 'deepvariant'
    germline: ['clair3','deepvariant']
    somatic: ['clairs','deepsomatic']
sv_calling:
    germline: ['sniffles']
    somatic: ['savana','severus','nanomonsv']

reference:
    build: GRCh38
    file: /path/to/reference.fasta

# Tool-specific annotations
annotation:
    blacklist: /path/to/data/delly/hg38-blacklist.v2.bed.gz
    pon: /path/to/data/severus/PoN_1000G_hg38.tsv.gz
    vntr: /path/to/data/severus/human_GRCh38_no_alt_analysis_set.trf.bed
    control_panel_path: /path/to/data/nanomonsv/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
    simple_repeat: /path/to/data/nanomonsv/simpleRepeat.bed.gz

params:
    trim: true
    read_qs: 10
    headcrop: 40
    tailcrop: 20
    minlen: 100
    minsvlen: 50

## Tools configuration
pepper:
  sif: /path/to/pepper.sif
  params: /path/to/data/ont_models/margin_phase
  
deepvariant:
  run_as_steps: true
  gpu: /path/to/deepvariant.gpu.sif
  cpu: /path/to/deepvariant.cpu.sif

deepsomatic: 
  gpu: /path/to/deepsomatic.gpu.sif
  cpu: /path/to/deepsomatic.cpu.sif

clair3:
  sif: /path/to/clair3.sif
  model: /path/to/data/ont_models/clair3_models/r1041_e82_400bps_sup_v520

clairS:
  sif: /path/to/clairS.sif
  platform: ont_r10_dorado_sup_5khz
  indel_calling: true
  subclone: true  # this enable subclone variant detection
```

## Usage

### Running the Complete Workflow
```bash
# Clone the repository
git clone https://github.com/bakeronit/nanopore_paired_tumour_workflow.git
cd nanopore_paired_tumour_workflow

# Configure workflow parameters
vi config/config.yaml
# Edit config.yaml with your paths and settings

# Check commands are valid using dry-run
snakemake --configfile config/config.yaml --directory /path/to/project/output --dry-run

# Run the workflow
nohup snakemake --configfile config/config.yaml --directory /path/to/project/output --profile $PWD/workflow/profiles/pbspro &
```

### Monitoring Progress (Optional)
```bash
# Check workflow status
snakemake --summary --directory /path/to/project/output --configfile config/config.yaml

# Generate workflow report
snakemake --report report.html --directory /path/to/project/output --configfile config/config.yaml
```

## Output Structure

Details of the output files can be found in [filegraph](docs/filegraph.png) plot

```bash
analysis/
├── bam/
│   ├── {sample}.bam
│   ├── {sample}.bam.bai
│   ├── {sample}.haplotagged.bam
│   └── {sample}.haplotagged.bam.bai
├── qc/
│   └── bam/
│       ├── {sample}.bamcov.txt
│       ├── {sample}.mosdepth.global.dist.txt
│       └── {sample}.mosdepth.summary.txt
├── snvs/
│   ├── clairs/
│   │  └── {sample_tumour}.{sample_normal}/
│   │       ├── indel.vcf.gz
│   │       ├── indel.vcf.gz.tbi
│   │       ├── output.vcf.gz
│   │       └── output.vcf.gz.tbi
│   │─ deepsomatic/
│   │   └── {sample_tumour}.{sample_normal}/
│   │       ├── output.vcf.gz
│   │       └── output.vcf.gz.tbi
│   │       ├── output.somatic.vcf.gz
│   │       └── output.somatic.vcf.gz.tbi
│   │── deepvariant/
│   │   └── {sample}/
│   │       ├── {sample}.vcf.gz
│   │       ├── {sample}.vcf.gz.tbi
│   │       ├── {sample}.g.vcf.gz
│   │       ├── {sample}.g.vcf.gz.tbi
│   │       ├── {sample}.passed.phased.vcf.gz
│   │       └── {sample}.passed.phased.vcf.gz.tbi
│   └── clair3/  #optional, may be excluded
├── svs/
│   ├── sniffles/
│   │   └── {sample}/
│   │       └── {sample}.vcf
│   ├── nanomonsv/
│   │   ├──{sample_tumour}.{sample_normal}/
│   │   │  ├── {sample_tumour}.{sample_normal}.nanomonsv.result.simple_repeat.svtype.passed.txt
│   │   │  ├── {sample_tumour}.{sample_normal}.nanomonsv.result.simple_repeat.svtype.txt
│   │   │  ├── {sample_tumour}.{sample_normal}.nanomonsv.result.txt
│   │   │  ├── {sample_tumour}.{sample_normal}.nanomonsv.result.vcf
│   │   │  ├── {sample_tumour}.{sample_normal}.nanomonsv.sbnd.result.txt
│   │   │  └──{sample_tumour}.{sample_normal}.nanomonsv.supporting_read.txt
│   │   └── {sample}/
│   │       ├── {sample}.bp_info.sorted.bed.gz
│   │       ├── {sample}.bp_info.sorted.bed.gz.tbi
│   │       ├── {sample}.deletion.sorted.bed.gz
│   │       ├── {sample}.deletion.sorted.bed.gz.tbi
│   │       ├── {sample}.insertion.sorted.bed.gz
│   │       ├── {sample}.insertion.sorted.bed.gz.tbi
│   │       ├── {sample}.rearrangement.sorted.bed.gz
│   │       └── {sample}.rearrangement.sorted.bed.gz.tbi
│   ├── delly/  #optional, may be excluded
│   ├── savana/
│   │   └── {sample_tumour}.{sample_normal}/
│   │       ├── {sample_tumour}.{sample_normal}.classified.somatic.bedpe
│   │       ├── {sample_tumour}.{sample_normal}.classified.somatic.vcf
│   │       ├── {sample_tumour}.{sample_normal}.classified.sv_breakpoints.vcf
│   │       ├── {sample_tumour}.{sample_normal}.inserted_sequences.fa
│   │       ├── {sample_tumour}.{sample_normal}.sv_breakpoints.bedpe
│   │       ├── {sample_tumour}.{sample_normal}.sv_breakpoints_read_support.tsv
│   │       └── {sample_tumour}.{sample_normal}.sv_breakpoints.vcf
│   └── severus/
│       └──{sample_tumour}.{sample_normal}/
│          ├── all_SVs/
│          ├── breakpoints_double.csv
│          ├── read_ids.csv
│          ├── read_qual.txt
│          ├── serverus_collapsed_dup.bed
│          ├── serverus.log
│          ├── serverus_LOG.bed
│          └── somatic_SVs/
│              ├── breakpoints_clusters_list.tsv
│              ├── breakpoints_clusters.tsv
│              ├── plots/
│              └── severus_somatic.vcf
├── cnvs/
│   └── savana/
│       └──{sample_tumour}.{sample_normal}/
│          ├── 10kbp_bin_ref_all_{sample_tumour}.{sample_normal}_with_SV_breakpoints.bed
│          ├── {sample_tumour}.{sample_normal}_allele_counts_hetSNPs.bed
│          ├── {sample_tumour}.{sample_normal}_fitted_purity_ploidy.tsv
│          ├── {sample_tumour}.{sample_normal}_ranked_solutions.tsv
│          ├── {sample_tumour}.{sample_normal}_raw_read_counts.tsv
│          ├── {sample_tumour}.{sample_normal}_read_counts_mnorm_log2r_segmented.tsv
│          └── {sample_tumour}.{sample_normal}_segmented_absolute_copy_number.tsv
├── mod/
│   ├── {sample}.bed.gz
│   ├── {sample}.bed.gz.tbi
│   ├── {sample}_1.bed.gz
│   ├── {sample}_1.bed.gz.tbi
│   ├── {sample}_2.bed.gz
│   ├── {sample}_2.bed.gz.tbi
│   ├── {sample}_ungrouped.bed.gz
│   ├── {sample}_ungrouped.bed.gz.tbi
│   └── {sample}.mod_summary.txt
└── reports/
    └── multiqc_report.html #TODO
```

### Key Output Files

- **Alignments**: `analysis/bam/{sample}.bam` - sorted genome alignment
- **Germline SNVs**: `analysis/snvs/deepvariant/{sample}/{sample}.passed.phased.vcf.gz`
- **Somatic SNVs**: `analysis/snvs/deepsomatic|clairS/{sample_tumour}.{sample_normal}/output.vcf.gz` - needs additional filters
- **Structural Variants**: `analysis/svs/{caller}/{sample}.vcf.gz`
- **QC Reports**: `analysis/qc/` or `analysis/qc/multiqc_report.html`
- **Methylation**: `analysis/mod/{sample}.bed.gz`

## Workflow Components

### Supported Tools

| Category | Tool | Purpose | 
|----------|------|---------|
| **Trimming** | [chopper](https://github.com/wdecoster/chopper) | Quality and adapter trimming |
| **Alignment** | [minimap2](https://github.com/lh3/minimap2) | Map reads to reference |
| **QC** | [mosdepth](https://github.com/brentp/mosdepth), [samtools/coverage](https://github.com/fbreitwieser/bamcov) | Coverage and quality metrics |
| **Germline SNVs** | [deepvariant](https://github.com/google/deepvariant), [whatshap](https://github.com/whatshap/whatshap) | Call germline variants and phasing |
| **Somatic SNVs** | [clairS](https://github.com/HKU-BAL/ClairS), [deepsomatic](https://github.com/google/deepsomatic) | Call somatic mutations |
| **Germline SVs** | [sniffles](https://github.com/fritzsedlazeck/Sniffles) | Call structural variants |
| **Somatic SVs** | [nanomonsv](https://github.com/friend1ws/nanomonsv), [savana](https://github.com/cortes-ciriano-lab/savana), [severus](https://github.com/KolmogorovLab/Severus) | Call somatic SVs |
| **Methylation** | [modkit](https://github.com/nanoporetech/modkit) | DNA methylation calling |
| **CNVs** | [savana](https://github.com/cortes-ciriano-lab/savana) | Call copy number alterations |
| **Merged SVs*** | [ssvanalyser](https://github.com/bakeronit/SSVAnalyser), [jasmine](https://github.com/mkirsche/Jasmine) | Merge structural variants |
| **Report*** | [multiqc](https://github.com/MultiQC/MultiQC) | Generate report |

*not yet implemented

### Analysis Modes

- **`germline`**: Process normal samples only for germline variant calling
- **`somatic`**: Paired tumour/normal analysis for somatic variants  
- **`all`**: Complete analysis including both germline and somatic variants

## Resource Requirements

### Compute Resources

In the `workflow/profiles/config.yaml`, you can change the resource requirements for each rule.
```yaml
set-threads:
  bc_summary: 1
  qc_bam_stats: 4
  qc_bam_cov: 1
  align_minimap2: 24
  merge_bam: 10
  ...
```

### Optimising the Workflow Resources Requsts on HPC

Although Snakemake can do benchmarking of resources used in each jobs in the workflow, the memory tracking is not the same as qstat. This workflow provides a script [qstat_jobs.py](workflow/scripts/qstat_jobs.py). At the end of the workflow, you could run `python workflow/scripts/qstat_jobs.py /path/to/work` to get a better idea of whether you are requesting too much or too little memory, CPUs, etc. This allows you to optimise your workflow resource requirements.

**An example of the `qstat_jobs.py` output suggests that I requested too much memory and too many CPUs:**
```bash
+-------------------+----------------+---------------------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| Job ID            | Job Name       | Wildcards                       | Host       |   CPU Req |   CPU Max |   CPU Avg |   Mem Req(GB) |   Mem Used(GB) | Time Req   | Time Used   |
+===================+================+=================================+============+===========+===========+===========+===============+================+============+=============+
| 30966588.hpcpbs02 | align_minimap2 | sample=COLO829_50, run=PAQ22155 | hpcnode068 |        24 |     13.21 |      0.24 |            20 |          18.85 | 02:00:00   | 00:18:23    |
+-------------------+----------------+---------------------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| 30966589.hpcpbs02 | align_minimap2 | sample=COLO829_50, run=PAQ45153 | hpcnode073 |        24 |      3.75 |      0.26 |            20 |          19.59 | 02:00:00   | 00:18:24    |
+-------------------+----------------+---------------------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| 30966590.hpcpbs02 | align_minimap2 | sample=COLO829_BL, run=PAQ23729 | hpcnode073 |        24 |      3.45 |      0.24 |            20 |          19.04 | 02:00:00   | 00:18:23    |
+-------------------+----------------+---------------------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| 30973416.hpcpbs02 | bamtohapmod    | sample=COLO829_BL               | hpcnode068 |         8 |      1.84 |      0.16 |            10 |           0.38 | 12:00:00   | 00:18:21    |
+-------------------+----------------+---------------------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
```

## Citation

If you use this workflow, please cite:

- **Snakemake**: Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012.
- **Individual tools**: Please cite the original publications for each tool used in the analysis

## License

This workflow is licensed under the MIT License. See [LICENSE](LICENSE) file for details.