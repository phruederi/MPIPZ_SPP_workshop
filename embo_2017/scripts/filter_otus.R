#!/projects/dep_psl/grp_psl/garridoo/tools/bin/Rscript

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# run script silently

sink(type="message")
options(warn=-1)

# parse arguments

args <- commandArgs(TRUE)
otu_table_file <- args[1]
output_file <- args[2]
threshold <- args[3]

# load normalized OTU table

otu_table <- read.table(args[1], sep="\t", header=T, check.names=F)
OTUId <- otu_table[, 1]
otu_table <- otu_table[, -1]

# thresholding

idx <- apply(otu_table * 100, 1, FUN=mean) >= threshold
otu_table_filtered <- otu_table[idx, ]
OTUId <- OTUId[idx]

# write filtered OTU table to the output file

otu_table_filtered <- data.frame(OTUId, otu_table_filtered)
write.table(otu_table_filtered, output_file, col.names=T, row.names=F, quote=F, sep="\t")

