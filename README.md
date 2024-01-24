## Snakemake workflow: nanopore\_paired\_tumour\_workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)

This workflow performs variant calling using nanopore data of tumour and paired normal samples.

### Tools:

**Basecalling** 
- [dorado](https://github.com/nanoporetech/dorado)
- [pod5](https://github.com/nanoporetech/pod5-file-format)

**Trimming**
- [chopper](https://github.com/wdecoster/chopper)

**Alignment**
- [minimap2](https://github.com/lh3/minimap2)

**QC**
- [mosdepth](https://github.com/brentp/mosdepth)
- [bamcov](https://github.com/fbreitwieser/bamcov)

**SNV calling**
- [clair3](https://github.com/HKU-BAL/Clair3)
- [clairS](https://github.com/HKU-BAL/ClairS)
- [deepvariant](https://github.com/google/deepvariant)
- [pepper-margin-deepvariant](https://github.com/kishwarshafin/pepper)
- [deepsomatic](https://github.com/google/deepsomatic)

**SV calling**
- [nanomonsv](https://github.com/friend1ws/nanomonsv)
- [delly](https://github.com/dellytools/delly)
- [savana](https://github.com/cortes-ciriano-lab/savana)
- [severus](https://github.com/KolmogorovLab/Severus)

**Mod call**
- [modkit](https://github.com/nanoporetech/modkit)

**Miscellaneous**
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
 
### Prerequisite 
**Genome**: reference assembly should be available and be indexed by minimap, samtools, and bwa.  

**Models**: model files for dorado, clair3, should be available and provided with path in `config/config.yaml`.  

**Annotation files**: Some tools require downloading extra data 
- nanomonsv
    - [Simple repeats](https://github.com/friend1ws/nanomonsv/wiki/An-example-on-removing-indels-within-simple-repeat).
    - [Control panel](https://zenodo.org/records/7017953)
    - [Misc scripts](https://github.com/friend1ws/nanomonsv/tree/master/misc)
- Severus
    - [VNTRs](https://github.com/KolmogorovLab/Severus/tree/main/vntrs)
 
### Sample

An example sample table (`config/samples.csv`):

|sample_id |run_id                         |flowcell_id|flowcell_version|donor_id|type  |raw_path                              |
|----------|-------------------------------|-----------|----------------|--------|------|--------------------------------------|
|COLO829_BL|Run folder name from seq centre|PAQXXXXX   |R10             |COLO829 |normal|/path/to/the/output/of/flowcell       |
|COLO829   |Run folder name from seq centre|PAQXXXXX   |R10             |COLO829 |tumour|/path/to/the/output/of/flowcell       |

*One flowcell a line. Sometimes a flowcell might be reloaded and generates two output, they should be in two lines.*