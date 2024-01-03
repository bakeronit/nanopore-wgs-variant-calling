## Snakemake workflow: nanopore\_paired\_tumour\_workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg)](https://snakemake.github.io)

This workflow performs variant calling using nanopore data of tumour and paired normal samples.


### Tools:

 **Basecalling**
 
 - [dorado](https://github.com/nanoporetech/dorado)

 **Alignment**

 - [minimap2](https://github.com/lh3/minimap2)

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


