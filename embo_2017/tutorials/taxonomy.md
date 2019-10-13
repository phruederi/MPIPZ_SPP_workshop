### Hands-on session - Taxonomic classification

#### Introduction

Classification is typically done only for OTU representative sequences instead of the full set of amplicon reads to save runtime. The representative sequences are produced by the OTU picking step in QIIME (morning session). Taxonomic classification must be combined with abundance information in order to generate estimates and plots of the taxonomic community structure. Abundance is the number of amplicon read within each OTU cluster.

#### Data

* `./data/tax/`: contains representative 16S sequences for classification
* `./data/tax/references/`: contains the taxonomy files and database 16S sequences used for classification
* `./data/tax/references/whole_genome/`: contains classifications generated from a genome assembly for comparison to 16S based assignments

#### Taxonomic classification of amplicon data

In QIIME, taxonomic classification is done using a wrapper script called [assign_taxonomy.py](http://qiime.org/scripts/assign_taxonomy.html), which can call different programs for classification. Here, we can make use of three different methods:

* a UCLUST-based search and assignment
* the RDP classifier implemented in mothur
* the SortMeRNA algorithm

For each method, one is free to choose a set of reference sequences and taxonomy. In this tutorial, you can select from

* `./data/tax/references/greengenes_13.8`: the [greengenes database](http://greengenes.secondgenome.com/) (97% identity clusters)
* `./data/tax/references/rdp_11.5`: [the Ribosomal Database Project](https://rdp.cme.msu.edu)
* `./data/tax/references/silva_128`: the [ARB-SILVA database](https://www.arb-silva.de)

#### Classifying amplicon reads

Pick the set of OTU representatives in `references/centroid.fna` and classify them by taxonomy using the QIIME methods. Refer to the [online documentation](http://qiime.org/scripts/assign_taxonomy.html) and specify an output folder for each combination of method and reference dataset.

#### Compare classifications

The output files are sorted differently. You can use the function `show_tax_assignments` in this tutorial to compare the results of different classification runs.

    show_tax_assignments output1.txt output2.txt | less -S

Try to spot the difference in the assignments. What does make a larger difference, a) the choice of classification method or b) the choice of reference data?

#### Online classification

Both ARB SILVA and RDP project have online classifiers on their website, which work much the same as the onces we run in this tutorial. The advantage is, that you don't need to install and run anything on the command line and have the latest database version readily available but you loose some control over the process and it can be more difficult to automate and reproduce the steps. There is usually also an upper file size limit, so it doesn't work for large samples.

* [RDP online classifier](https://rdp.cme.msu.edu/seqmatch/seqmatch_intro.jsp)
* [SINA online alignment and classifiation](https://www.arb-silva.de/aligner/)

You can trim your input FASTA files using the `seqkit` command. For instance, to extract the first 10 16S seqences, do

    seqkit head -n 10 < original.fna > extracted.fna

To select a particular sequences, run

    seqkit grep name_of_16s_seq < original.fna > extracted.fna

You might need to transfer these files to your local computer before uploading.

Explore and compare the results from the offline and online classification.

#### Continue from the morning session

Try to continue from the morning session with the OTU clusters that you generated and classify the cluster representatives. See the documentation for `make_otu_table.py` and `summarize_taxa.py` in QIIME to ammend your OTU tables with taxonomic annotation for further analysis like diversity estimation or plotting.

#### Whole-genome classifications

Metagenome and draft genome assemblies contain more information for taxonomic classification than just 16S sequences. In particular, strains may deviate in genome sequence but be similar in 16S sequences. Therefore, it is worth comparing them as well. `references/whole_gnome/wgs.tax` contains whole-genome classifications in the same format as the 16S results which were generated using the program [taxator-tk](https://github.com/fungs/taxator-tk). Compare to the 16S assignments using the function `how_tax_assignments` as above.

Genome assemblies are not as easy to work with as single gene sequences, because they may be very fragmented and may be contaminated with sequences from other genomes. For instance, the genome assembly of the leaf326 OTU (`references/whole_genome/leaf326.fna`) consists of 108 contigs with a total size of 4.6 Mb. When this assembly is subjected to classification, some contigs will be informative while others will yield very general taxon predictions. The output of a taxator-tk classification for this assembly can be viewed by running

    less -S references/whole_genome/leaf326.taxpath.tsv

A consensus taxonomic prediction can be generated from these results as well as hypothetical contamination estimated.

