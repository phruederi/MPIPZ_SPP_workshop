### Hands-on session - Alpha and beta diversity estimates

Different sampling depths (varying number of reads per sample) can have an
important effect on the calculation of diversity estimates. For instance,
very deep samples will appear to correspond to more complex communities
than samples with a smaller number of sequences. There are two main ways
in which this problem is addressed: rarefaction and normalization.

There is some controversy in the field as to what approach is preferable,
(see e.g. McMurdie *et al.*, 2014 and Weiss *et al.*, 2017). The answer
very likely depends on particular factors of each experiment and caution
should be exercised before deciding which technique should be applied.

#### Rarefaction

One way to artificially make all samples of equal depth consists on randomly
picking an equal number of reads from each sample and discarding the rest. This
number can be given by the smallest sample size or a higher number of reads,
in which case samples under that threshold will be discarded entirely. This
approach is highly encouraged for alpha diversity estimates but can lead
to complications for other diversity analyses or analyses of differentially
abundant OTUs, in addition to the fact that generally a considerable percentage
of the reads are thrown away.

This step can be easily performed using QIIME by running the following
command:

``single_rarefaction.py -i otu_table.biom -o otu_table_raref.biom -d 500``

for an equal sampling depth of 500 reads per sample.

For compatibility, the rarefied OTU table in biom format can be converted
to a text file by typing:

``biom convert -i otu_table_raref.biom --table-type="OTU table" --to-tsv -o otu_table_raref.txt``

#### Normalization (total sum)

One way to normalize OTU tables (per samples) while avoiding throwing away
useful data is to divide each count of each OTU by the total sum of all reads
found in a given sample, which is known as total sum normalization. This can
be done using a simple R script as follows:

``./scripts/normalize_otu_table.R otu_table.txt otu_table_norm.txt``

``biom convert -i otu_table_norm.txt -o otu_table_norm.biom --table-type="OTU table" --to-json``

#### Filtering of low abundance OTUs

One additional step which might be particularly relevant in case of large
differences in sample size is to filter very low abundant OTUs. The reason
for this is two fold: first, since a large proportion of low abundance OTUs
are artifacts that have not been successfully filtered, this step ensures
further removal of errors to the expense of a (small) loss of depth. Secondly,
this step can help to partially address the problem of inflated alpha-diversity
in deeper samples with only a small trade-off.

There are various ways in which this step can be done, for example, by filtering
all OTUs which are not present above a certain threshold in at least one treatment
or condition (for highly varying conditions) or by removing all OTUs whose
mean relative abundance in the entire dataset falls bellow the threshold (for
studies with more homogeneous conditions). In general, filtering OTUs with
mean abundances bellow a reasonably low value (e.g. 0.1% or 0.001%) is advisable.
This can be done using an R script by typing:

``./scripts/filter_otus.R otu_table_norm.txt otu_table_norm_filtered.txt 0.01``

``biom convert -i otu_table_norm_filtered.txt -o otu_table_norm_filtered.biom --table-type="OTU table" --to-json``

Finally, we can remove some comments from the text format OTU tables to
ensure compatibility with scripts used in downstream analyses

``sed -i '/# Const/d;s/#OTU ID.//g;s/OTUId.//g' otu_table*.txt``

#### Calculating indices of alpha-diversity (within sample diversity)

Measures of alpha diversity attempt to provide a quantitative estimate of the
'complexity' of a (microbial) community. There are various ways in which
complexity might play a role, for example, one sample can have a very large
number of species but be dominated by one in particular, whereas another
can contain a smaller number of species but at roughly equal abundances.
This can lead to varying estimates, depending on what features we take into
account. It is therefore a good idea to calculate and check different
indices (e.g. Shannon, Chao, number of species).
It is easy to calculate these within QIIME by applying the following commands
to the rarefied OTU table:

``alpha_diversity.py -i otu_table_raref.biom -o shannon.txt -m shannon``

``alpha_diversity.py -i otu_table_raref.biom -o chao.txt -m chao1``

``alpha_diversity.py -i otu_table_raref.biom -o observed_otus.txt -m observed_species``

#### Calculating measures of beta-diversity (between sample diversity)

Beta diversity metrics provide a quantitative evaluation of how much samples
differ from one another. They provide a framework for extracting features
from microbial community data which would not be apparent by looking at
the distributions of species relative abundances of each sample. In many
cases, it is possible to perform statistical analyses on these metrics
as we shall see later on. There are multiple distances that can be roughly
grouped into quantitative (take into account OTU abundances), for example
Bray Curtis or UniFrac or qualitative (they consider only presence / absence)
such as Jaccard or unweighted UniFrac.

It is easy to obtain matrices of pairwise distances using QIIME and
taking as input always normalized (or rarefied) OTU tables. For instance
to obtain Bray Curtis distance matrices, we can use the following commands:

``beta_diversity.py -i otu_table_norm.biom -m bray_curtis -o ./``

``mv bray_curtis_otu_table_norm.txt bray_curtis.txt``

``sed -i 's/^\t//g' bray_curtis.txt``

By to Bray Curtis, UniFrac distances take into account the phylogeny (relatedness)
of the OTUs found in the dataset in order to quantify beta diversity. The rational
behind this is the assumption that two samples that differ in the presence of
two very closely related species are more similar than two samples that vary in
the presence of very diverse OTUs. In order to estimate UniFrac distances,
therefore, we need first to infer a phylogeny or species tree for the OTU
representative sequences.

First, we have to align the FASTA file containing the representatives. This can
be done using a variety of methods, for example ClustalO:

``clustalo --seqtype=DNA --threads 1 -i rep_seqs.fasta -o rep_seqs_aligned.fasta``

And subsequently use QIIME script to generate the tree from the multiple
sequence alignment, by running the command:

``make_phylogeny.py -i rep_seqs_aligned.fasta -o rep_seqs.tree``

Finally, we can use the same script as above to calculate UniFrac distances,
both weighted

``beta_diversity.py -i otu_table_norm.biom -m weighted_unifrac -o ./ -t rep_seqs.tree``

``mv weighted_unifrac_otu_table_norm.txt weighted_unifrac.txt``

``sed -i 's/^\t//g' weighted_unifrac.txt``

and unweighted

``beta_diversity.py -i otu_table_norm.biom -m unweighted_unifrac -o ./ -t rep_seqs.tree``

``mv unweighted_unifrac_otu_table_norm.txt unweighted_unifrac.txt``

``sed -i 's/^\t//g' unweighted_unifrac.txt``

These matrices of pairwise distances can be used in downstream analyses,
for example to perform dimensionality reduction and plotting or examine the
cluster structure of the dataset.

