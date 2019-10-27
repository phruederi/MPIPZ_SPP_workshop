##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFlg22 <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 3, each = 6)) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}

save_pheatmap_jpg <- function(x, filename, width = 480, height = 480) {
  jpeg(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width = 8, height = 8) {
  cairo_pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('magrittr')
suppressPackageStartupMessages(library('tximport'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('RColorBrewer'))
library('pheatmap')

outdir <- '~/RNA-Seq_visual'

if(!file.exists(outdir)) {
  dir.create(outdir)
} else {}

anno <- read_csv('/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_index/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
mutate(Description = Description %>% {if_else(is.na(.), '', .)})

##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step1: load alignments')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

wd <- '/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_align'
setwd(wd)

annoSample <- tibble(condition = paste0(rep(c('Col0', 'flg22', 'MeJA'), each = 6),
                                        '_s_',
                                        rep(c('rep1_', 'rep2_', 'rep3_'), 6),
                                        rep(c('batch1', 'batch2'), each = 3) %>% rep(3))) %>%
  mutate(sample = paste0(condition, '_ath_kallisto'))

kres <- file.path(wd, annoSample$sample, 'abundance.h5') %>%
  set_names(annoSample$condition) %>%
  tximport(type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step2: normalization counts')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

condi <- c('Col0', 'flg22', 'MeJA')
sampleTable <- data.frame(condition = factor(rep(condi, each = 6), levels = condi))
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
## Q1: how to transform row counts to log2-scaled form? ##
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~PCA plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step3: PCA plot')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(3, name = 'Dark2'))

## Q2: which function is used for PCA ##
pca <- rld %>%
  assay %>%
  t %>%
  prcomp

percentVar <- pca$sdev^2 %>%
  {./sum(.)} %>%
  {100 * .} %>%
  round
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, condition = colData(rld)[, 1], ID = rownames(colData(rld)), Batch = rep(c('batch1', 'batch2'), each = 3) %>% rep(3))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3, alpha = 0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols)) +
  geom_text(aes(label = Batch), hjust = -0.3, vjust = 0) +
  xlim(c(-40, 40)) +
  ylim(c(-40, 40))
ggsave(plot = p, file.path(outdir, 'PCA_singleend_batch.jpg'))
ggsave(plot = p, file.path(outdir, 'PCA_singleend_batch.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

