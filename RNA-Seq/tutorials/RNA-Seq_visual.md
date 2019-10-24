## SPP workshop - Expression data visualization in R
###  October 2019 | MPIPZ Cologne

<!-- content start -->

- [1. Prerequisites](#1-prerequisites)
- [2. Visualization](#2-visualization)
    - [2.1 PCA plot](#21-pca-plot)
    - [2.2 Heatmap](#22-heatmap)
    - [2.3 Run visualization](#23-run-visualization)
- [3. Practice](#3-practice)
- [References](#references)
    
<!-- content end -->

## 1. Prerequisites

1. Alignment results

We will use the whole pair-end alignments by kallisto for demonstration and single-end alignments by kallisto for practice. 

> pair-end alignments by kallisto

```bash
├── Col0_p_rep1_ath_kallisto
├── Col0_p_rep2_ath_kallisto
├── Col0_p_rep3_ath_kallisto
├── flg22_p_rep1_ath_kallisto
├── flg22_p_rep2_ath_kallisto
├── flg22_p_rep3_ath_kallisto
```

> single-end alignments by kallisto

```bash
├── Col0_s_rep1_ath_kallisto
├── Col0_s_rep2_ath_kallisto
├── Col0_s_rep3_ath_kallisto
├── flg22_s_rep1_ath_kallisto
├── flg22_s_rep2_ath_kallisto
├── flg22_s_rep3_ath_kallisto
├── MeJA_s_rep1_ath_kallisto
├── MeJA_s_rep2_ath_kallisto
└── MeJA_s_rep3_ath_kallisto
```

2. Scripts

Copy the scripts we need for visualization.

```
cp /netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_scripts/workshop_visual* ~
```

## 2. Visualization

### 2.1 PCA plot 

### 2.2 Heatmap

### 2.3 Run visualization

## 3. Practice

```bash
/netscratch/common/MPIPZ_SPP_workshop/software/R-3.6.1/bin/Rscript \
~/workshop_visual_pairend.R
```

## References

* Love MI, Huber W, Anders S. **Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.** *Genome Biol.* 2014;15(12):550.

* Castrillo G\*, Teixeira PJ\*, Paredes SH\*, Law TF, de Lorenzo L, Feltcher ME, Finkel OM, Breakfield NW, Mieczkowski P, Jones CD, Paz-Ares J, Dangl JL. **Root microbiota drive direct integration of phosphate stress and immunity.** *Nature* 2017;543(7646):513-518. Single-end data RNA-Seq data.

* Völz R\*, Kim SK, Mi J, Rawat AA, Veluchamy A, Mariappan KG, Rayapuram N, Daviere JM, Achard P, Blilou I, Al-Babili S, Benhamed M, Hirt H\*. **INDETERMINATE-DOMAIN 4 (IDD4) coordinates immune responses with plant-growth in Arabidopsis thaliana.** PLoS Pathog. 2019;15(1):e1007499. Pair-end RNA-Seq data.

