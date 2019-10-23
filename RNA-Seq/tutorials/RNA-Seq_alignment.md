## SPP workshop - Read alignment and pseudo-alignment
###  October 2019 | MPIPZ Cologne

<!-- content start -->

- [1. Prerequisites](#1-prerequisites)
- [2. Alignment by HISAT2](#2-alignment-by-hisat2)
    - [2.1 Index](#21-index)
    - [2.2 Alignment](#21-alignment)
- [3. Pseudo-alignment by kallisto](#3-pseudo-alignment-by-kallisto)
    - [3.1 Index](#31-index)
    - [3.2 Alignment and quantification](#32-alignment-and-quantification)
- [References](#references)
    
<!-- content end -->

## 1. Prerequisites

1. Clean fastq files

We will use pair-end data for demonstration and single-end data for practice.

> pair-end clean data has two conditions: Col0 and Col0 treated with flg22. Each condition has three replicates.

```
├── Col0_p_rep1_small_R1.fq.gz
├── Col0_p_rep1_small_R2.fq.gz
├── Col0_p_rep2_small_R1.fq.gz
├── Col0_p_rep2_small_R2.fq.gz
├── Col0_p_rep3_small_R1.fq.gz
├── Col0_p_rep3_small_R2.fq.gz
├── flg22_p_rep1_small_R1.fq.gz
├── flg22_p_rep1_small_R2.fq.gz
├── flg22_p_rep2_small_R1.fq.gz
├── flg22_p_rep2_small_R2.fq.gz
├── flg22_p_rep3_small_R1.fq.gz
├── flg22_p_rep3_small_R2.fq.gz
```

> single-end clean data set has three conditions: Col0, Col0 treated with flg22, and Col0 treated with methyl jasmonate (MeJA). Each condition has three replicates.

```
├── Col0_s_rep1_small.fq.gz
├── Col0_s_rep2_small.fq.gz
├── Col0_s_rep3_small.fq.gz
├── flg22_s_rep1_small.fq.gz
├── flg22_s_rep2_small.fq.gz
├── flg22_s_rep3_small.fq.gz
├── MeJA_s_rep1_small.fq.gz
├── MeJA_s_rep2_small.fq.gz
├── MeJA_s_rep3_small.fq.gz
```

2. Scripts.

Copy the scripts we need for alignment.

```bash
cp /netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_scripts/workshop_align* ~
```

## 2. Alignment by HISAT2

### 2.1 Index

[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) makes an index of ***Arabidopsis* genomes [TAIR10](ftp://ftp.ensemblgenomes.org/pub/plants/release-45/fasta/arabidopsis_thaliana/dna/)(`Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`). The index can be used for multiple RNA-Seq experiments, as long as the organism is the same. Building an HISAT2 index is time-consuming, and we have already built one looks like:

```bash
RNA-Seq_index/athht2index/
├── genome.1.ht2
├── genome.2.ht2
├── genome.3.ht2
├── genome.4.ht2
├── genome.5.ht2
├── genome.6.ht2
├── genome.7.ht2
└── genome.8.ht2
```

The script for building the index is:

```bash
## do not run the codes below！！
mkdir athht2index
cd athht2index
${HISAT2_PATH}/hisat2-build -f \
              ${REF_PATH}/ath/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
              genome
```

### 2.2 Alignment



## 



## References

* Kim D, Langmead B, Salzberg SL. **HISAT: a fast spliced aligner with low memory requirements.** *Nat Methods.* 2015;12(4):357-60.

* Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter. **Near-optimal probabilistic RNA-seq quantification.** *Nature Biotechnology* 2016;34(5):525–27.

* Sahraeian SME, Mohiyuddin M, Sebra R, Tilgner H, Afshar PT, Au KF, Bani Asadi N, Gerstein MB, Wong WH, Snyder MP, Schadt E, Lam HYK. **Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis.** *Nat Commun.* 2017;8(1):59.

* 
