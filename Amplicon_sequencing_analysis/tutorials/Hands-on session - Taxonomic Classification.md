# Hands-on session - Taxonomic Classification


**Introduction**
Classification is typically done only for OTU representative sequences/ASVs instead of the full set of amplicon reads to save runtime. Taxonomic classification must be combined with abundance information in order to generate estimates and plots of the taxonomic community structure.

**Data**
- ./reference_data/silva_132/silva_132_97_16S.fna: contains representative 16S sequences for classification
- ./reference_data/silva_132/taxonomy_7_levels.txt: contains the taxonomy files and database 16S sequences used for classification
- asv-seqs.qza: contains the representative ASV sequences to be classified
If you perform OTU clustering and get OTU representative sequences, you can use the following command to convert the fasta format to qza format:
```
qiime tools import --input-path otu_result/otus.fa --output-path otu_result/otus-seqs.qza --type FeatureData[Sequence]
```

**Taxonomic classification of amplicon data**
In QIIME2, taxonomic classification is done using a plugin feature-classifier, which can call different programs for classification. Here, we can make use of three different methods:
- BLAST+: local sequence alignment followed by consensus taxonomy classification
*Local alignment*
```
5' ACTACTAGATTACTTACGGATCAGGTACTTTAGAGGCTTGCAACCA 3' 
             |||| |||||| |||||||||||||||
          5' TACTCACGGATGAGGTACTTTAGAGGC 3'
```
- VSEARCH: global sequence alignment followed by consensus taxonomy classification
*Global alignment*
```
5' ACTACTAGATTACTTACGGATCAGGTACTTTAGAGGCTTGCAACCA 3'
   |||||||||||    |||||||  |||||||||||||| |||||||
5' ACTACTAGATT----ACGGATC--GTACTTTAGAGGCTAGCAACCA 3'
```
- classify-sklearn: Naive Bayes Classifier

For each method, one is free to choose a set of reference sequences and taxonomy. In general, there are three databases to choose, including [Greengenes], [RDP] and [Silva]. But only Silva keeps updating. In this tutorial, we only introduce how to classify reads taxonomically based on Silva database.
- ./data/tax/silva_132: the ARB-SILVA database (Download here: https://ftp.arb-silva.de/qiime/)


**Classifying amplicon reads**
Pick the set of OTU/ASV representatives in references/centroid.fna and classify them by taxonomy using the QIIME2 methods. Refer to the online documentation and specify an output folder for each combination of method and reference dataset.
----

Import the reference dataset into qiime2:
*If you download the silva database from the above-mentioned website, you need to modify the taxonomy file (taxonomy_7_levels.txt) by adding "Feature ID[tab]Taxon" to the first line.*
```
qiime tools import --input-path ./data/tax/silva_132/silva_132_97_16S.fna --output-path ./data/tax/silva_132/silva_132_97_16S.qza --type FeatureData[Sequence]
qiime tools import --input-path ./data/tax/silva_132/taxonomy_7_levels.txt --output-path ./data/tax/silva_132/taxonomy_7_levels.qza --type FeatureData[Taxonomy]
```

Run Blast+ to classify reads taxonomically
```
mkdir taxonomy
qiime feature-classifier classify-consensus-blast --i-query asv_result/asv-seqs.qza --i-reference-reads ./reference_data/silva_132/silva_132_97_16S.qza --i-reference-taxonomy ./reference_data/silva_132/taxonomy_7_levels.qza --o-classification taxonomy/silva_blast_taxonomy.qza
qiime tools export --input-path taxonomy/silva_blast_taxonomy.qza --output-path taxonomy/silva_blast_taxonomy
```
Now you can go into the *silva_blast_taxonomy* folder to check the annotation file:
```
cd taxonomy/silva_blast_taxonomy
head taxonomy.tsv
cd ../..
```

----


Run Vsearch to classify reads taxonomically
```
qiime feature-classifier classify-consensus-vsearch --i-query asv_result/asv-seqs.qza --i-reference-reads ./reference_data/silva_132/silva_132_97_16S.qza --i-reference-taxonomy ./reference_data/silva_132/taxonomy_7_levels.qza --o-classification taxonomy/silva_vsearch_taxonomy.qza
qiime tools export --input-path taxonomy/silva_vsearch_taxonomy.qza --output-path taxonomy/silva_vsearch_taxonomy
```
Also check the output:
```
head taxonomy/silva_vsearch_taxonomy/taxonomy.tsv
```

----


Run Bayes to classify reads taxonomically. Firstly you need to train the reference data to get the likehood probabilities.

Train the classifier (Calculate the likelihood of the Bayesian formula)
```
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ./reference_data/silva_132/silva_132_97_16S.qza --i-reference-taxonomy ./reference_data/silva_132/taxonomy_7_levels.qza --o-classifier ./reference_data/silva_132/silva_132_ref_97_nb_classifier.qza
```
Taxonomic assignments of reads
```
qiime feature-classifier classify-sklearn --i-reads asv_result/asv-seqs.qza --i-classifier ./reference_data/silva_132/silva_132_ref_97_nb_classifier.qza --o-classification taxonomy/silva_nb_taxonomy.qza
qiime tools export --input-path taxonomy/silva_nb_taxonomy.qza --output-path taxonomy/silva_nb_taxonomy
```
----
**Visualize community structure**
You can plot barplots to visualize the taxonomic community structure at different taxonomic levels. In order to reach this, you have to collapse the data at specified taxonomic level first. For example, if you want plot the community structure at family level, collapsing means summming up the abundance of OTUs/ASVs annotated to the same families. In QIIME2, there is a plugin designated to collapse the profiling table at different taxonomic level, but each time, only one taxonomic level can be specified. You can collapse the data at family by running the commands:
```
qiime taxa collapse --i-table asv_result/asv_table.qza --i-taxonomy taxonomy/silva_nb_taxonomy.qza --p-level 5 --o-collapsed-table asv_result/asv_table_family.qza
qiime tools export --input-path asv_result/asv_table_family.qza --output-path asv_result/collapsed_table
mv asv_result/collapsed_table/feature-table.biom asv_result/collapsed_table/asv_table_family.biom
biom convert -i asv_result/collapsed_table/asv_table_family.biom --to-tsv -o asv_result/collapsed_table/asv_table_family.tsv
```
2 - phylum | 3 - class |  4 - order | 5 - family | 6 - genus
Check the file:
```
head asv_result/collapsed_table/asv_table_family.tsv
```
Run the command to visualize the community structure:
```
Rscript ./script/plot_community_structure.R asv_result/collapsed_table/asv_table_family.tsv family ./asv_result/collapsed_table/
```
This command will generate 3 pdf files. One plots all the families observed in the data; the other two files filter out the families whose relative abundances are smaller than 0.1% in all samples and whose averagee relative abundances are smaller than 1% in all samples respetively.




[Greengenes]: <https://greengenes.secondgenome.com/>
[RDP]: <https://rdp.cme.msu.edu/>
[Silva]: <https://www.arb-silva.de>