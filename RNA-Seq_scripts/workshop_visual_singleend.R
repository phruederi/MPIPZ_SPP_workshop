##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
    res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
    return(res)
}

checkFlg22 <- function(v, threshold) {

    require('magrittr')

    res <- v %>%
    split(rep(1 : 3, each = 3)) %>%
    sapply(checkZeros, threshold) %>%
    all

    return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('magrittr')
suppressPackageStartupMessages(library('tximport'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('RColorBrewer'))

outdir <- '~/RNA-Seq_visual'

if(!file.exists(outdir)) {
  dir.create(outdir)
} else {}

anno <- read_csv('/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_index/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
mutate(Description = Description %>% {if_else(is.na(.), '', .)})

##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_align_data'
setwd(wd)

annoSample <- tibble(condition = paste0(rep(c('Col0', 'flg22', 'MeJA'), each = 3),
                                     '_s_',
                                     rep(c('rep1', 'rep2', 'rep3'), 3))) %>%
  mutate(sample = paste0(condition, '_ath_kallisto'))

kres <- file.path(wd, annoSample$sample, 'abundance.h5') %>%
  set_names(annoSample$condition) %>%
  tximport(type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
condi <- c('Col0', 'flg22', 'MeJA')
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|x and |0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ] %>%
  DESeq

## count transformation
rld <- rlog(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~PCA plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(3, name = 'Dark2'))

pca <- rld %>%
  assay %>%
  t %>%
  prcomp
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols))

ggsave(plot = p, file.path(outdir, 'PCA_singleend.jpg'))
ggsave(plot = p, file.path(outdir, 'PCA_singleend.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

