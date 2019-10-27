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
    split(rep(1 : 2, each = 3)) %>%
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

wd <- '/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_align_data'
setwd(wd)

annoSample <- tibble(condition = paste0(rep(c('Col0', 'flg22'), each = 3),
                                     '_p_',
                                     rep(c('rep1', 'rep2', 'rep3'), 2))) %>%
  mutate(sample = paste0(condition, '_ath_kallisto'))

kres <- file.path(wd, annoSample$sample, 'abundance.h5') %>%
  set_names(annoSample$condition) %>%
  tximport(type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step2: normalization counts')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

condi <- c('Col0', 'flg22')
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
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~PCA plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step3: PCA plot')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(3, name = 'Dark2')[1:2])

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
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, condition = colData(rld)[, 1], ID = rownames(colData(rld)))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols))

ggsave(plot = p, file.path(outdir, 'PCA_pairend.jpg'))
ggsave(plot = p, file.path(outdir, 'PCA_pairend.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step3: PCA plot')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
cond <- list(c('Col0', 'flg22'))

## DEGs
resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs =  list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Col0_p_rep1 : Col0_vs_flg22_log2FoldChange) %>%
  arrange(Col0_vs_flg22_padj)

## selected DEGs
heatmapData <- degres %>%
  rownames %>%
  match(res$ID[1:50], .) %>%
  assay(ntd)[., ]

df <- colData(degres) %>% as.data.frame
heatmapPlot <- pheatmap(heatmapData,
                        cluster_rows = TRUE,
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        annotation_col = df,
                        annotation_colors = list(condition = c(Col0 = '#1B9E77', flg22 = '#D95F02')))
save_pheatmap_jpg(heatmapPlot, file.path(outdir, 'heatmap_pairend.jpg'))
save_pheatmap_pdf(heatmapPlot, file.path(outdir, 'heatmap_pairend.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
