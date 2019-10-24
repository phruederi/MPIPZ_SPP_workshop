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

annoSample <- tibble(condition = paste0(rep(c('Col0', 'flg22', 'MeJA'), each = 3),
                                     '_s_',
                                     rep(c('rep1', 'rep2', 'rep3'), 3))) %>%
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
cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(3, name = 'Dark2'))

pca <- rld %>%
  assay %>%
  t %>%
  prcomp
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, condition = colData(rld)[, 1], ID = rownames(colData(rld)))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols))

ggsave(plot = p, file.path(outdir, 'PCA_singleend.jpg'))
ggsave(plot = p, file.path(outdir, 'PCA_singleend.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step3: PCA plot')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
cond <- list(c('Col0', 'flg22'),
             c('Col0', 'MeJA'))

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
  select(ID, Gene : Description, Col0_s_rep1 : Col0_vs_MeJA_log2FoldChange) %>%
  arrange(Col0_vs_flg22_padj)

## selected DEGs
sDEGs <- c(res$ID[1:30],
           res %>% arrange(Col0_vs_MeJA_padj) %>% .$ID %>% .[1:30]) %>%
  unique

heatmapData <- degres %>%
  rownames %>%
  match(sDEGs, .) %>%
  assay(ntd)[., ]

df <- colData(degres) %>% as.data.frame
heatmapPlot <- pheatmap(heatmapData,
                        cluster_rows = TRUE,
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        annotation_col = df,
                        annotation_colors = list(condition = c(Col0 = '#1B9E77', flg22 = '#D95F02', MeJA = '#7570B3')))
save_pheatmap_jpg(heatmapPlot, file.path(outdir, 'heatmap_singleend.jpg'))
save_pheatmap_pdf(heatmapPlot, file.path(outdir, 'heatmap_singleend.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


