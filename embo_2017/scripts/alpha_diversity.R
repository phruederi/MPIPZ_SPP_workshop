
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")

# directories

results.dir <- "../data/diversity/"
figures.dir <- "../figures/"

# files

design.file <- paste(results.dir, "design.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_norm.txt", sep="")

shannon.file <- paste(results.dir, "shannon.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=F, fill=T)

shannon <- read.table(shannon.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[idx, idx]

# remove non-bacterial OTUs

taxonomy <- taxonomy[taxonomy[, 2]=="Bacteria", ]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove negative control samples

neg_control_samples <- design$SampleID[design$compartment=="negative_control"]
idx <- !colnames(otu_table) %in% neg_control_samples
otu_table <- otu_table[, idx]
design <- design[idx, ]

# remove individual nodules samples

ind_nod_samples <- design$SampleID[design$compartment=="individual_nodule"]
idx <- !colnames(otu_table) %in% ind_nod_samples
otu_table <- otu_table[, idx]
design <- design[idx, ]
 
# subset samples of interest from distance matrices

idx <- shannon[, 1] %in% design$SampleID

shannon <- shannon[idx, ]

### alpha diversity

colors <- data.frame(group=c("pooled_nodules", "root", "rhizosphere", "soil"),
                     color=c(c_red, c_very_dark_green, c_dark_red, c_dark_brown))

shapes <- data.frame(group=c("wiltype", "mutant", "soil"),
                     shape=c(19, 0, 18))

# shannon index

index <- cbind(shannon[, 2], design[match(shannon[, 1], design$SampleID), ])
colnames(index)[1] <- "value"

index$genotype <- factor(index$genotype, levels=shapes$group)
index$compartment <- factor(index$compartment, levels=colors$group)

# reorder boxplots

l <- c("soil", "rhizosphere", "root", "pooled_nodules")
index$compartment <- factor(index$compartment, levels=l)
colors <- colors[match(l, colors$group), ]

p <- ggplot(index, aes(x=compartment, y=value, color=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
            geom_jitter(aes(shape=genotype), position=position_jitter(0.17), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="", y="shannon index") +
            main_theme

ggsave(paste(figures.dir, "shannon.pdf", sep=""), p)

