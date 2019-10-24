## SPP workshop - Read alignment and pseudo-alignment
###  October 2019 | MPIPZ Cologne

<!-- content start -->

- [1. Prerequisites](#1-prerequisites)
- [2. Index](#2-index)
- [3. Alignment](#3-alignment)
    - [3.1 Alignment by HISAT2](#31-alignment-by-hisat2)
    - [3.2 Pseudo-alignment by kallisto](#32-pseudo-alignment-by-kallisto)
    - [3.3 Run (pseudo-)alignment](#33-run-pseudo-alignment)
- [4. Practice](#4-practice)
- [References](#references)
    
<!-- content end -->

## 1. Prerequisites

1. Clean fastq files

We generated a set of small fastq files (randomly select ~1k reads) to accelerate the analysis process. We will use pair-end data for demonstration and single-end data for practice.

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

> <span id="single-end-data">single-end</span> clean data set has three conditions: Col0, Col0 treated with flg22, and Col0 treated with methyl jasmonate (MeJA). Each condition has three replicates.

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

Copy the scripts we need for alignment through `ssh`.

```bash
cp /netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_scripts/workshop_align* ~
```

## 2. Index

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

The script for building the index looks like:

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

The following scripts shows how to align `Col0_p_rep1_small_R1.fq.gz` and `Col0_p_rep1_small_R2.fq.gz`

```bash
## do not run the codes below！！
## -p: threads number.
## -x: HISAT2 index path.
## -1 and -1: paired fastq files.
## -S: output path.
hisat2 -p 2 \
       -x athht2index/genome \
       -1 Col0_p_rep1_small_R1.fq.gz \
       -2 Col0_p_rep1_small_R2.fq.gz \
       -S Col0_p_rep1_small_ath_hisat2.sam
```

### 3.2 Pseudo-alignment by kallisto

```bash
## do not run the codes below！！
## -t: threads number.
## -i: kallisto index.
## -o: output path.
kallisto quant \
         -t 2 \
         -i ath.kindex \
         -o Col0_p_rep1_small_ath_kallisto \
         Col0_p_rep1_small_R1.fq.gz \
         Col0_p_rep1_small_R2.fq.gz 
```

### 3.3 Run (pseudo-)alignment

1. Run pair-end (pseudo-)alignment through `ssh`:

```bash
## it may take 5 min
bash ~/workshop_align_pairend.sh
```

2. Check alignment information

The output looks like:

```
====================================
Kallisto using ath cDNA for Col0_p_rep1_small.

[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 48,359
[index] number of k-mers: 45,567,218
[index] number of equivalence classes: 93,973
[quant] running in paired-end mode
[quant] will process pair 1: /netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_smallclean_data/Col0_p_rep1_small_R1.fq.gz
                             /netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_smallclean_data/Col0_p_rep1_small_R2.fq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 967 reads, 949 reads pseudoaligned
[quant] estimated average fragment length: 197.992
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 477 rounds

HISAT2 using ath genome for Col0_p_rep1_small.
967 reads; of these:
  967 (100.00%) were paired; of these:
    59 (6.10%) aligned concordantly 0 times
    890 (92.04%) aligned concordantly exactly 1 time
    18 (1.86%) aligned concordantly >1 times
    ----
    59 pairs aligned concordantly 0 times; of these:
      1 (1.69%) aligned discordantly 1 time
    ----
    58 pairs aligned 0 times concordantly or discordantly; of these:
      116 mates make up the pairs; of these:
        75 (64.66%) aligned 0 times
        41 (35.34%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
96.12% overall alignment rate
[bam_sort_core] merging from 0 files and 2 in-memory blocks...
=====================================
```

## 4. Practice

1. Open and fix the script `workshop_align_singleend_practice.sh` at your home folder (enter home folder by execute `cd ~`). Three places marked with `Q1 - Q3` needed to be fixed. 

> `Q1`: pattern for [single-end fastq data](#single-end-data)?

> `Q2`: command for kallisto alignment?

> `Q3`: command for samtools sort reads?

2. Test your codes through `ssh`:

```
## 1. run single-end (pseudo-)alignment, it may take 5 min
bash ~/workshop_align_singleend_practice.sh

## 2. overview of alignment results
tree ~/tree RNA-Seq_align_data
```

A fixed script `workshop_align_singleend_answer.sh` at your home folder can be used as a reference.

## References

* Kim D, Langmead B, Salzberg SL. **HISAT: a fast spliced aligner with low memory requirements.** *Nat Methods.* 2015;12(4):357-60.

* Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter. **Near-optimal probabilistic RNA-seq quantification.** *Nature Biotechnology* 2016;34(5):525–27.

* Sahraeian SME, Mohiyuddin M, Sebra R, Tilgner H, Afshar PT, Au KF, Bani Asadi N, Gerstein MB, Wong WH, Snyder MP, Schadt E, Lam HYK. **Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis.** *Nat Commun.* 2017;8(1):59.

* Castrillo G\*, Teixeira PJ\*, Paredes SH\*, Law TF, de Lorenzo L, Feltcher ME, Finkel OM, Breakfield NW, Mieczkowski P, Jones CD, Paz-Ares J, Dangl JL. **Root microbiota drive direct integration of phosphate stress and immunity.** *Nature* 2017;543(7646):513-518. Single-end data RNA-Seq data.

* Völz R\*, Kim SK, Mi J, Rawat AA, Veluchamy A, Mariappan KG, Rayapuram N, Daviere JM, Achard P, Blilou I, Al-Babili S, Benhamed M, Hirt H\*. **INDETERMINATE-DOMAIN 4 (IDD4) coordinates immune responses with plant-growth in Arabidopsis thaliana.** PLoS Pathog. 2019;15(1):e1007499. Pair-end RNA-Seq data.
