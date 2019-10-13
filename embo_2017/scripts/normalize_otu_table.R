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

# load count OTU table

otu_table <- read.table(otu_table_file, sep="\t", header=T, check.names=F)
OTUId <- otu_table[, 1]
otu_table <- otu_table[, -1]

# total sum normalization

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

# write normalized OTU table to the output file

otu_table_norm <- data.frame(OTUId, otu_table_norm)
write.table(otu_table_norm, output_file, col.names=T, row.names=F, quote=F, sep="\t")

