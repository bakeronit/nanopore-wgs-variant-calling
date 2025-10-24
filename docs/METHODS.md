# Long Read Sequencing Variant Calling

**The Medical Genomics Team@QIMRB**

> This document is generated with **Claude AI** (without curation) to provide an overview of the methods used in the pipeline.

## Sample Processing and Sequencing

Paired tumour and normal tissue samples were sequenced using Oxford Nanopore Technologies whole-genome sequencing on PromethION flow cells (R10.4.1 chemistry). Basecalling was performed upstream using the Dorado basecaller with the super-accurate model.

## Bioinformatics Analysis Pipeline

All bioinformatic analyses were performed using a Snakemake workflow (version ≥8.20.1) with containerized tools executed via Singularity/Apptainer on a high-performance computing cluster. The complete workflow is available at https://github.com/MedicalGenomicsLab/nanopore_paired_tumour_workflow.

### Read Preprocessing and Quality Control

Raw sequencing reads in unaligned BAM (uBAM) format underwent quality filtering and trimming using chopper. Specifically, the first 40 bp and last 20 bp of each read were trimmed (headcrop=40, tailcrop=20), and reads shorter than 100 bp or with mean quality score below Q10 were removed. Following trimming, base modification tags (MM tags) were repaired using modkit repair to ensure consistency with the trimmed read coordinates.

### Read Alignment

Quality-filtered reads were aligned to the GRCh38 human reference genome (GCA_000001405.15_GRCh38_no_alt_analysis_set) using minimap2 (v2.24) with the lr:hq preset for high-quality long reads. The alignment was performed with MD tag generation (-Y --MD flags) and phasing information retention (-y flag). Aligned reads were sorted using samtools (v1.17) and indexed. For samples sequenced across multiple flow cells, BAM files were merged using samtools merge.

### Alignment Quality Assessment

Alignment quality metrics were generated using mosdepth (v0.3.3) to calculate depth distribution and coverage statistics across the genome, and samtools coverage to compute per-chromosome coverage metrics.

### Germline Variant Calling

#### Small Variants (SNVs and Indels)

Germline small variants were called using two complementary approaches:

1. **DeepVariant** (v1.9.0): Variant calling was performed in three stages using the ONT R10.4 model:
   - **make_examples**: Generated pileup images from aligned reads with alt-aligned pileup using diff_channels mode, phased read sorting, and read quality filtering (minimum MAPQ=5). A small variant model was employed to improve indel detection with custom thresholds (indel GQ≥17, SNP GQ≥9).
   - **call_variants**: Applied the trained deep learning model to classify variants from the generated examples.
   - **postprocess_variants**: Generated final VCF and gVCF files with variant quality scores.

2. **Clair3** (v1.1.2): Variants were called using the r1041_e82_400bps_sup_v520 model with phasing enabled and including all contigs (standard and non-standard chromosomes).

#### Variant Phasing

High-quality germline variants (FILTER=PASS) from DeepVariant were phased using WhatsHap (v1.7), which leverages read-backed phasing to assign variants to haplotype blocks. The phasing was performed using the aligned BAM files with the --ignore-read-groups flag to pool all reads together.

#### Haplotagging

Aligned reads were assigned to haplotypes using WhatsHap haplotag based on the phased VCF files. This process added HP (haplotype) tags to reads, assigning them to haplotype 1 or 2 based on the variant phasing information. Haplotagged BAM files were generated for both tumour and normal samples and indexed for downstream analyses.

### Somatic Variant Calling

#### Small Variants (SNVs and Indels)

Somatic small variants were identified using two paired tumour-normal callers:

1. **DeepSomatic** (v1.9.0): Somatic variants were called chromosome-by-chromosome using the ONT model to manage computational resources efficiently. The analysis processed chromosomes 1-22, X, Y, and M separately with 24 shards per chromosome. Resulting VCF files were concatenated using bcftools concat, filtered for homozygous variants (GT="1/1") and PASS quality, then sorted to generate the final somatic variant set.

2. **ClairS** (v0.4.3): Paired somatic variant calling was performed using the ont_r10_dorado_sup_5khz platform setting with the r1041_e82_400bps_sup_v520 Clair3 model. Indel calling was enabled (--enable_indel_calling), and subclonal variant detection was activated (--enable_verdict) to identify variants present at lower variant allele frequencies.

### Structural Variant Calling

#### Germline Structural Variants

Germline structural variants (SVs) were detected using Sniffles (v2.0) with a minimum SV length threshold of 50 bp. Sniffles identifies SVs by analyzing split-read and read-depth signatures from long-read alignments.

#### Somatic Structural Variants

Somatic SVs were identified using three complementary algorithms to maximize detection sensitivity and specificity:

1. **nanomonsv** (v0.8.0): This tool operates in two stages:
   - **parse**: Extracts potential breakpoints, deletions, insertions, and rearrangements from individual BAM files, generating sorted BED/BEDPE files for each sample.
   - **get**: Compares tumour and matched normal breakpoint signatures to identify somatic SVs. The analysis included single breakend detection (--single_bnd), local assembly refinement using racon (--use_racon), small indel detection (minimum 10 bp), and filtering against a panel of normals from the Human Pangenome Reference Consortium (HPRC Year 1 nanopore data, n=30 samples). Results were post-processed to annotate simple repeat regions and classify insertion sequences.

2. **SAVANA** (v1.3.0): Somatic SV calling was performed in two steps:
   - **run**: Identified candidate SV breakpoints from paired tumour-normal BAM files with minimum SV length of 50 bp, single breakend detection enabled, and a minimum of 3 supporting reads required.
   - **classify**: Applied a machine learning classifier to distinguish somatic from germline SVs, with copy number alteration (CNA) rescue mode enabled (--cna_rescue) and ONT-specific parameters (--ont).

3. **Severus** (v1.1): Detected somatic SVs using haplotagged BAM files and phased VCF files from both tumour and normal samples. The analysis included variable number tandem repeat (VNTR) annotation, minimum SV size of 50 bp, single breakpoint detection (--single-bp), between-junction insertion detection (--between-junction-ins), and loss of heterozygosity (LOH) identification (--output-LOH). Collapsed duplications and supporting read IDs were also output for further investigation.

### Copy Number Alteration Analysis

Somatic copy number alterations (CNAs) were inferred using SAVANA cna, which integrates multiple data sources:
- Read depth signals from tumour and normal samples binned at 10 kbp resolution
- Allelic imbalance from heterozygous SNPs identified in the phased germline VCF
- Structural variant breakpoints to inform copy number transitions
- ENCODE blacklist regions (hg38 v2) to exclude problematic genomic regions

The tool performed segmentation of normalized log2 ratios, estimated tumour purity and ploidy, and calculated absolute copy numbers for genomic segments. Multiple solutions were ranked based on fit quality.

### DNA Methylation Analysis

DNA methylation patterns were profiled from the base modification tags present in the aligned BAM files using modkit (v0.5.0):

1. **Combined methylation calling**: CpG methylation was quantified genome-wide using modkit pileup with strand combination (--combine-strands) and CpG-specific analysis (--cpg), generating bedMethyl format files compressed and indexed with bgzip and tabix.

2. **Haplotype-resolved methylation**: Methylation calls were partitioned by haplotype using the haplotagged BAM files with the --partition-tag HP option, generating separate methylation profiles for haplotype 1 and haplotype 2, plus ungrouped reads.

3. **Modification summary**: Global modification statistics were computed using modkit summary to quantify the proportion of modified bases across all read contexts.

### Workflow Management and Reproducibility

The entire analysis workflow was orchestrated using Snakemake (v8.20.1), enabling parallel job execution, automatic dependency resolution, and resource management on a PBS Pro high-performance computing cluster. All software tools were containerized using Singularity/Apptainer to ensure reproducibility across different computing environments. Computational resource usage was monitored and optimized using custom scripts that parsed PBS job statistics to identify over- or under-provisioned jobs.

### Reference Data and Annotations

- **Reference genome**: GRCh38 (GCA_000001405.15_GRCh38_no_alt_analysis_set)
- **Panel of normals (Severus)**: 1000 Genomes Project samples for hg38
- **Panel of normals (nanomonsv)**: HPRC Year 1 nanopore sequencing data (n=30)
- **Blacklist regions**: ENCODE Blacklist v2 for hg38
- **Tandem repeats**: TRF-based VNTR annotations for GRCh38
- **Simple repeats**: Pre-built simple repeat annotations for nanomonsv

### Data Availability

The workflow implementation, configuration files, and containerized software environments are publicly available at https://github.com/MedicalGenomicsLab/nanopore_paired_tumour_workflow.
