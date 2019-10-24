# Hands-on Session - Alpha and beta diversity estimates

Different sampling depths (varying number of reads per sample) can have an important effect on the calculation of diversity estimates. For instance, very deep samples will appear to correspond to more complex communities than samples with a smaller number of sequences. There are two main ways in which this problem is addressed: rarefaction and normalization.

There is some controversy in the field as to what approach is preferable, (see e.g. McMurdie et al., 2014 and Weiss et al., 2017). The answer very likely depends on particular factors of each experiment and caution should be exercised before deciding which technique should be applied.

**Rarefaction**
One way to artificially make all samples of equal depth consists on randomly picking an equal number of reads from each sample and discarding the rest. This number can be given by the smallest sample size or a higher number of reads, in which case samples under that threshold will be discarded entirely. This approach is highly encouraged for alpha diversity estimates but can lead to complications for other diversity analyses or analyses of differentially abundant OTUs, in addition to the fact that generally a considerable percentage of the reads are thrown away.

This step can be easily performed using QIIME2 by running the following command:
```
qiime feature-table rarefy --i-table asv_result/asv_table.qza --p-sampling-depth 500 --o-rarefied-table asv_result/asv_table_raref.qza
```
for an equal sampling depth of 500 reads per sample.

**Normalization(total sum)**
One way to normalize OTU tables (per samples) while avoiding throwing away useful data is to divide each count of each OTU by the total sum of all reads found in a given sample, which is known as total sum normalization. This can be done using a simple R script as follows:
```
./scripts/normalize_otu_table.R asv_result/asv_table_rename.tsv asv_result/asv_table_norm.tsv
```

**Filtering of low abundance OTUs**
One additional step which might be particularly relevant in case of large differences in sample size is to filter very low abundant OTUs. The reason for this is two fold: first, since a large proportion of low abundance OTUs are artifacts that have not been successfully filtered, this step ensures further removal of errors to the expense of a (small) loss of depth. Secondly, this step can help to partially address the problem of inflated alpha-diversity in deeper samples with only a small trade-off.

There are various ways in which this step can be done, for example, by filtering all OTUs which are not present above a certain threshold in at least one treatment or condition (for highly varying conditions) or by removing all OTUs whose mean relative abundance in the entire dataset falls bellow the threshold (for studies with more homogeneous conditions). In general, filtering OTUs with mean abundances bellow a reasonably low value (e.g. 0.1% or 0.001%) is advisable. This can be done using an R script by typing:
```
./scripts/filter_otus.R asv_result/asv_table_norm.tsv asv_result/asv_table_norm_filtered.tsv 0.01
```
Then you can import the filtered asv table into QIIME2 again by running the following two commands (Note: Only biom file in hdf5 format is compatible with *qiime tools import*): 
```
biom convert -i asv_result/asv_table_norm_filtered.tsv -o asv_result/asv_table_norm_filtered.biom --table-type="OTU table" --to-hdf5
qiime tools import --input-path asv_result/asv_table_norm_filtered.biom --output-path asv_result/asv_table_norm_filtered.qza --type FeatureTable[Frequency]
```


**Calculating indices of alpha-diversity (within sample diversity)**
Measures of alpha diversity attempt to provide a quantitative estimate of the 'complexity' of a (microbial) community. There are various ways in which complexity might play a role, for example, one sample can have a very large number of species but be dominated by one in particular, whereas another can contain a smaller number of species but at roughly equal abundances. This can lead to varying estimates, depending on what features we take into account. It is therefore a good idea to calculate and check different indices (e.g. Shannon, Chao, number of species). It is easy to calculate these within QIIME2 by applying the following commands to the rarefied OTU table:
```
mkdir alpha_diversity
qiime diversity alpha --i-table asv_result/asv_table_raref.qza --p-metric shannon --o-alpha-diversity alpha_diversity/shannon.qza
qiime diversity alpha --i-table asv_result/asv_table_raref.qza --p-metric chao1 --o-alpha-diversity alpha_diversity/chao1.qza
qiime diversity alpha --i-table asv_result/asv_table_raref.qza --p-metric osd --o-alpha-diversity alpha_diversity/observed.qza
qiime tools export --input-path alpha_diversity/shannon.qza --output-path alpha_diversity
```

**Calculating measures of beta-diversity (between sample diversity)**
Beta diversity metrics provide a quantitative evaluation of how much samples differ from one another. They provide a framework for extracting features from microbial community data which would not be apparent by looking at the distributions of species relative abundances of each sample. In many cases, it is possible to perform statistical analyses on these metrics as we shall see later on. There are multiple distances that can be roughly grouped into quantitative (weighted / take into account OTU abundances), for example Bray Curtis or UniFrac or qualitative (unweighted / they consider only presence / absence) such as Jaccard or unweighted UniFrac.

It is easy to obtain matrices of pairwise distances using QIIME2 and taking as input always rarefied (or normalized) OTU tables. For instance to obtain Bray Curtis distance matrices, we can use the following commands:
```
mkdir beta_diversity
qiime diversity beta --i-table asv_result/asv_table_raref.qza --p-metric braycurtis --p-n-jobs 1 --o-distance-matrix beta_diversity/bc_distance.qza
qiime tools export --input-path beta_diversity/bc_distance.qza --output-path beta_diversity
```

Different to Bray Curtis, UniFrac distances take into account the phylogeny (relatedness) of the OTUs found in the dataset in order to quantify beta diversity. The rational behind this is the assumption that two samples that differ in the presence of two very closely related species are more similar than two samples that vary in the presence of very diverse OTUs. In order to estimate UniFrac distances, therefore, we need first to infer a phylogeny or species tree for the OTU representative sequences.

First, we have to generate a phylogenetic tree for the representatives by running the command:
```
mkdir phylo_tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences asv_result/asv-seqs.qza --p-n-threads 1 --o-alignment phylo_tree/alignment.qza --o-masked-alignment phylo_tree/masked_alignment.qza --o-tree phylo_tree/unrooted_tree.qza --o-rooted-tree phylo_tree/rooted_tree.qza
```
Then we can calculate weighted Unifrac distance with the command:
```
qiime diversity beta-phylogenetic --i-table asv_result/asv_table_raref.qza --i-phylogeny phylo_tree/rooted_tree.qza --p-metric weighted_unifrac --p-n-jobs 1 --o-distance-matrix beta_diversity/weighted_unifrac_distance.qza
qiime tools export --input-path beta_diversity/weighted_unifrac_distance.qza --output-path beta_diversity
```
You can also run unweighted Unifrac with the command:
```
qiime diversity beta-phylogenetic --i-table asv_result/asv_table_raref.qza --i-phylogeny phylo_tree/rooted_tree.qza --p-metric unweighted_unifrac --p-n-jobs 1 --o-distance-matrix beta_diversity/nweighted_unifrac_distance.qza
```