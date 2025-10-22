# Nanopore whole-genome sequencing variant calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.20.5-brightgreen.svg)](https://snakemake.github.io) [![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3100/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Sylabs SIF](https://img.shields.io/badge/Container-Sylabs%20SIF-5E72E4?style=flat-square)](https://cloud.sylabs.io/library/jiazhang/workflows/long_read_wgs_pipeline)

A Snakemake workflow for germline and somatic variant calling from Oxford Nanopore long-read whole-genome sequencing data of paired tumour/normal samples.

## Overview

This workflow processes Oxford Nanopore long-read WGS sequencing data to identify:
- **Germline variants**: SNVs and SVs present in normal tissue
- **Somatic variants**: SNVs and SVs acquired in tumour tissue
- **Copy number alterations**: Using multiple SV callers
- **Phased variants**: With haplotagging support

## Installation and Setup

### Requirements
- **Snakemake ≥8.20.1**
- **Singularity/Apptainer** (for containerized tools)
- **Python 3.11+**

### Container Setup

#### Option 1: Use Pre-built Container (Recommended)
```bash
# Download the main workflow container
apptainer pull library://jiazhang/workflows/long_read_wgs_pipeline:1.0
```

#### Option 2: Build Container Locally if new versions are needed
```bash
# Build the container from definition file
apptainer build long_read_wgs_pipeline.sif long_read_wgs.def
```

### Download Required Models and Annotations

#### Reference Genome

A constant reference genome and indexed genome is required for all samples.

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
- [CliarS](https://github.com/HKU-BAL/ClairS) `--platform` parameter changes with ONT models, need to keep it up-to-date. 

This four tools are not included in the pipeline container, but each tool has a containerized version available:

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

 - Variable number tandem repeats: [VNTRs](https://github.com/KolmogorovLab/Severus/tree/main/vntrs) annotations for common references, customed annotations can be created using [findTandemRepeats](https://github.com/PacificBiosciences/pbsv/blob/master/annotations/findTandemRepeats)
 - [PON](https://github.com/KolmogorovLab/Severus/tree/main/pon), used for **tumour only** cases. 


## Configuration

### Sample Sheets

Create `config/samples.csv` with your sample information:

```csv
sample_id,flowcell_id,flowcell_version,donor_id,type
sample_t,PAQ23729,R10,PATIENT01,tumour
sample_n,PAQ24563,R10,PATIENT01,normal
```

> for QIMR Promethion sequencing, this sample files can be generated as the [nanopore basecalling pipeline](https://github.com/bakeronit/nanopore_dorado_basecalling), and manually add the type column.

**Required columns:**
- `sample_id`: Unique identifier for each sample
- `flowcell_id`: Flowcell barcode (e.g., PAQ23729)
- `flowcell_version`: Flowcell chemistry (R9/R10, best supporting R10)
- `donor_id`: Patient/donor identifier (must match between tumour/normal pairs)
- `type`: Sample type (`normal` or `tumour`, if running in germline mode, then samples are processed as simple sample, no somatic calling)

### Workflow Configuration

Edit `config/config.yaml`:

```yaml
# Sample and analysis configuration
samples: /path/to/samples.csv
run_mode: "all"  # Options: "germline", "somatic", "all"
pull_containers: false  # Set true if containers need pulling, better not do this on HPC

# Reference genome
reference:
  build: GRCh38
  file: /path/to/reference.fasta

# Variant calling tools to run
snv_calling:
  phased_snv_from: 'deepvariant'
  germline: ['deepvariant']
  somatic: ['clairs','deepsomatic']

sv_calling:
  germline: ['sniffles']
  somatic: ['savana','severus','nanomonsv']

# Container paths (update with your paths)
clair3:
  sif: /path/to/clair3_v1.1.2.sif
  model: /path/to/clair3_models/r1041_e82_400bps_sup_v520

deepvariant:
  cpu: /path/to/deepvariant_1.9.0.sif
  gpu: /path/to/deepvariant_v1.9.0-gpu.sif

# Additional annotations
annotation:
  blacklist: /path/to/hg38-blacklist.v2.bed.gz
  pon: /path/to/PoN_1000G_hg38.tsv.gz
  vntr: /path/to/human_GRCh38_no_alt_analysis_set.trf.bed
  control_panel_path: /path/to/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
  simple_repeat: /path/to/simpleRepeat.bed.gz
```

## Usage

### Running the Complete Workflow

```bash
# Clone the repository
git clone https://github.com/bakeronit/nanopore_paired_tumour_workflow.git
cd nanopore_paired_tumour_workflow

# Configure workflow parameters
cp config/config.yaml.template config/config.yaml
# Edit config.yaml with your paths and settings

# Check commands are valid using dry-run
snakemake --configfile config/config.yaml --directory /path/to/project/output --dry-run

# Run the workflow
nohup snakemake --configfile config/config.yaml --directory /path/to/project/output --profile workflow/profiles/pbspro &
```

### Monitoring Progress (Optional)

```bash
# Check workflow status
snakemake --summary --directory /path/to/project/output --configfile config/config.yaml

# Generate workflow report
snakemake --report report.html --directory /path/to/project/output --configfile config/config.yaml
```

## Output Structure

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

- **Alignments**: `analysis/bam/{sample}.bam`  # sorted genome alignment
- **Germline SNVs**: `analysis/snvs/deepvariant/{sample}/{sample}.passed.phased.vcf.gz`
- **Somatic SNVs**: `analysis/snvs/deepsomatic|clairS/{sample_tumour}.{sample_normal}/output.vcf.gz`  # Need to have more filters
- **Structural Variants**: `analysis/svs/{caller}/{sample}.vcf.gz`
- **QC Reports**: `analysis/qc/` or `analysis/qc/multiqc_report.html`
- **Methylation**: `analysis/mod/{sample}.bed.gz`

## Workflow Components

### Supported Tools

| Category | Tool | Purpose | 
|----------|------|---------|
| **Trimming** | chopper | Quality and adapter trimming |
| **Alignment** | minimap2 | Map reads to reference |
| **QC** | mosdepth, samtools/coverage | Coverage and quality metrics |
| **Germline SNVs** | deepvariant, whatshap | Call germline variants, and phasing |
| **Somatic SNVs** | clairS, deepsomatic | Call somatic mutations |
| **Germline SVs** | sniffles | Call structural variants |
| **Somatic SVs** | nanomonsv, savana, severus | Call somatic SVs |
| **Methylation** | modkit | DNA methylation calling |
| **CNVs** | savana | Call copy number alterations |
| **Merged SVs*** | SSVanalyser, Jasmine | Merge structural variants |
| **Report*** | multiqc | Generate report |

**not yet implemented*

### Analysis Modes

- **`germline`**: Process normal samples only for germline variant calling
- **`somatic`**: Paired tumour/normal analysis for somatic variants  
- **`all`**: Complete analysis including both germline and somatic variants

## Resource Requirements

### Compute Resources

In the `workflow/profiles/config.yaml`, you can change the resources requirements for each rule.

```yaml
set-threads:
  bc_summary: 1
  qc_bam_stats: 4
  qc_bam_cov: 1
  align_minimap2: 24
  merge_bam: 10
  ...
```

### Keep optimising the workflow resources requirements
Although snakemake can do benchmarking of the workflow, the memory tracking is not the same as qstat. To get better estiamte, at the end of the workflow, you could run `python workflow/scripts/qstat_jobs.py work/` to get a better idea of whether you are requsting too much or little memory, CPUs, etc. So you can optimise your workflow resources requirements.

**An example of the `qstat_jobs.py` output indicating that I requsted too much memory and too many CPUs)**

```bash
+-------------------+-------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| Job ID            | Job Name          | Host       |   CPU Req |   CPU Max |   CPU Avg |   Mem Req(GB) |   Mem Used(GB) | Time Req   | Time Used   |
+===================+===================+============+===========+===========+===========+===============+================+============+=============+
| 30952074.hpcpbs02 | coral_reconstruct | hpcnode062 |        16 |      1.22 |      0.81 |            52 |          18.09 | 48:00:00   | 01:30:43    |
+-------------------+-------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
| 30953787.hpcpbs02 | coral_reconstruct | hpcnode068 |        16 |      0.99 |      0.96 |            52 |          12.37 | 48:00:00   | 00:42:05    |
+-------------------+-------------------+------------+-----------+-----------+-----------+---------------+----------------+------------+-------------+
```

## Citation

If you use this workflow, please cite:

- **Snakemake**: Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012.
- **Individual tools**: Please cite the original publications for each tool used in the analysis

## License

This workflow is licensed under the MIT License. See `LICENSE` file for details.