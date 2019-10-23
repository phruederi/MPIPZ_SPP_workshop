## SPP workshop - Read alignment and pseudo-alignment
###  October 2019 | MPIPZ Cologne

<!-- content start -->

- [1. Prerequisite](#1-prerequisite)
- [2. Alignment by HISAT2](#2-alignment-by-hisat2)
    - [2.1 Index](#21-index)
    - [2.2 Alignment](#21-alignment)
- [3. Pseudo-alignment by kallisto](#3-pseudo-alignment-by-kallisto)
    - [3.1 Index](#31-index)
    - [3.2 Alignment and quantification](#32-alignment-and-quantification)
- [References](#references)
    
<!-- content end -->

## 1. Prerequisite

We need to prepare two things before alignment:

1. clean fastq files

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

2. 



## References

* Kim D, Langmead B, Salzberg SL. **HISAT: a fast spliced aligner with low memory requirements.** *Nat Methods.* 2015;12(4):357-60.

* Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter. **Near-optimal probabilistic RNA-seq quantification.** *Nature Biotechnology* 2016;34(5):525–27.

* Sahraeian SME, Mohiyuddin M, Sebra R, Tilgner H, Afshar PT, Au KF, Bani Asadi N, Gerstein MB, Wong WH, Snyder MP, Schadt E, Lam HYK. **Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis.** *Nat Commun.* 2017;8(1):59.

* 
