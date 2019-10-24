## SPP workshop - Read alignment and pseudo-alignment
###  October 2019 | MPIPZ Cologne

<!-- content start -->

- [1. Prerequisites](#1-prerequisites)
- [2. Index](#2-index)
- [3. Alignment](#3-alignment)
    - [3.1 Alignment by HISAT2](#31-alignment-by-hisat2)
    - [3.2 Pseudo-alignment by kallisto](#32-pseudo-alignment-by-kallisto)
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

## 2. Index

### 2.1 Index

[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) at first makes an index according to the *Arabidopsis* [TAIR10](http://plants.ensembl.org/Arabidopsis_thaliana/Info/Index) genome (`Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`), whereas [kallisto](https://pachterlab.github.io/kallisto/) uses the transcriptome (`Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz`). The index can be used for multiple RNA-Seq experiments, as long as the organism is the same. Building index is time-consuming, and we have already built one looks like:

```bash
RNA-Seq_index
├── athht2index
│   ├── genome.1.ht2
│   ├── genome.2.ht2
│   ├── genome.3.ht2
│   ├── genome.4.ht2
│   ├── genome.5.ht2
│   ├── genome.6.ht2
│   ├── genome.7.ht2
│   └── genome.8.ht2
├── ath.kindex
```

The script for building the index:

```bash
## do not run the codes below！！
## HISAT2 index
mkdir athht2index
cd athht2index
hisat2-build -f \
  Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  genome

## kallisto index
kallisto index \
  -i ath.kindex \
  Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
```

## 3. Alignment

### 3.1 Alignment by HISAT2

### 3.2 Pseudo-alignment by kallisto

## References

* Kim D, Langmead B, Salzberg SL. **HISAT: a fast spliced aligner with low memory requirements.** *Nat Methods.* 2015;12(4):357-60.

* Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter. **Near-optimal probabilistic RNA-seq quantification.** *Nature Biotechnology* 2016;34(5):525–27.

* Sahraeian SME, Mohiyuddin M, Sebra R, Tilgner H, Afshar PT, Au KF, Bani Asadi N, Gerstein MB, Wong WH, Snyder MP, Schadt E, Lam HYK. **Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis.** *Nat Commun.* 2017;8(1):59.

* Castrillo G\*, Teixeira PJ\*, Paredes SH\*, Law TF, de Lorenzo L, Feltcher ME, Finkel OM, Breakfield NW, Mieczkowski P, Jones CD, Paz-Ares J, Dangl JL. **Root microbiota drive direct integration of phosphate stress and immunity.** *Nature* 2017;543(7646):513-518. Single-end data RNA-Seq data.

* Völz R\*, Kim SK, Mi J, Rawat AA, Veluchamy A, Mariappan KG, Rayapuram N, Daviere JM, Achard P, Blilou I, Al-Babili S, Benhamed M, Hirt H\*. **INDETERMINATE-DOMAIN 4 (IDD4) coordinates immune responses with plant-growth in Arabidopsis thaliana.** PLoS Pathog. 2019;15(1):e1007499. Pair-end RNA-Seq data.
