# Hands-on session - Pre-processing of raw amplicon data
----
**Introduction**

This pipeline is used to preprocess amplicon sequencing data from natural or synthetic communities. Both QIIME2 and Usearch are used here. [QIIME2] can be downloaded here https://docs.qiime2.org/2019.7/install/. [Usearch] has two versions: 32X and 64X. 32X is freely available, but it's limited to low memory usage. 64X is charged without any memory limitation. For natural communities, steps including demutiplexing, reads merge and OTU clustering or ASV denoising are introduced. As for synthetic communities, we mainly introduce how to use exact matches to profile the data.
Here we mainly introduced how to use QIIME2 to preprocess the data.

----

**File Description**
In general, there are 4 files needed for analyze amplicon data as follows:
- *Forward.fastq.gz*: The sequencing reads from the forward strand
- *Reverse.fastq.gz*: The sequencing reads from the reverse strand
- *Barcode.fastq.gz*: The barcode reads at each coordination
- *Mapping.txt*: A metadata file includes a detailed description of the samples.

In this case, only one barcode is attached to each sample, but sometimes double barcodes are attacthed to a sample, which can increase the number of pooled samples during sequencing library construction.

----
# Natural Communities
----
### Pipeline for OTU clustering
**Activate QIIME2**
To make sure that the python compiler finds all the necessary packages, we need to activate the QIIME2 installation environments.
```
conda activate qiime2-2019.7
```

**Step1: Import sequencing data into QIIME2**
QIIME2 has a specific file format (.qza) to install information, so we need to convert the data into a QIIME readable format. Here we only show how to import data in which only one barcode is attached to each sample. An example how to import data with double barcodes is given in the extended material.
```
mkdir raw_data 
mv *.gz raw_data
cd raw_data
mv Barcode.fastq.gz barcodes.fastq.gz
mv Forward.fastq.gz forward.fastq.gz
mv Reverse.fastq.gz reverse.fastq.gz
cd ..
```
Here we create a folder that exclusively deposits the three sequencing data files and rename the three sequencing files as QIIME2 can only recognize the files with those names.

```
qiime tools import --type EMPPairedEndSequences --input-path ./raw_data --output-path ./data.qza
```

**Step2: Sample Demutiplexing**
Using the information that matches sample IDs with their corresponding barcodes (file Mapping.txt), it is possible to assign each read to a sample (and discard those erroneous barcodes)
```
qiime demux emp-paired --m-barcodes-file mapping.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --i-seqs data.qza --o-per-sample-sequences data_demux.qza --o-error-correction-details data_demux-details.qza
```
Stastistic of data in each sample
```
qiime demux summarize --i-data data_demux.qza --o-visualization data_demux.qzv
```
Convert .qza file format into a human-readable format. Now you can open the website https://view.qiime2.org/, drag the file into the browser and visualize the statistical summary of each sample.
```
qiime tools export --input-path data_demux.qza --output-path demultiplexed_data/
```
Then we can deactivate the QIIME2 environment at this moment:
```
conda deactivate
```


**Step3: Merge Pair-end Reads and Perform Quality Filtering**
This step is intended to remove low-quality reads in each sample
```
awk '{if(NR>1){print $0}}' mapping.txt| cut -f 1 | while read line;do F1=`ls demultiplexed_data/*L001_R1_001.fastq.gz | grep ${line}_`; R1=`ls demultiplexed_data/*L001_R2_001.fastq.gz | grep ${line}_`; flash2 $F1 $R1 -M 250 -o $line -x 0.25 -d demultiplexed_data/;done
```
You can check the output from flash2:
```
ls demultiplexed_data/
```
Remove the intermediate files:
```
rm demultiplexed_data/*.gz
rm demultiplexed_data/*notCombined_*.fastq
rm demultiplexed_data/*.hist*
```
Quality filtering
```
for i in demultiplexed_data/*.extendedFrags.fastq; do name=`basename $i .extendedFrags.fastq`; ../software/usearch/usearch -fastq_filter $i -fastaout demultiplexed_data/${name}.trimmed.fasta -fastq_maxee 3 -fastq_maxns 0 -relabel $name.;done
```
Rename the reads in each sample in order to make the names compatible in the script generating the OTU table.
```
for i in demultiplexed_data/*.trimmed.fasta;do name=`basename $i .trimmed.fasta`; awk -v sample="$name" 'BEGIN{c=1}{if(/^>/){print ">"sample"."c";barcodelabel="sample;c++} else{print $0}}' $i > demultiplexed_data/$name.trimmed_renamed.fasta;done
rm demultiplexed_data/*.trimmed.fasta
```


**Step4: Dereplicate Reads**
Sequences without errors that come from the same biological entity are identical. Ideally, most of the dataset would consist of many replicates of the same real template sequences. Although this is unfortunately not the case (~50% of the reads only, due to errors and artifacts), it is possible to substantially speed-up downstream processing by removing replicates (while keeping count), e.g. by using the following command in USEARCH:
```
cat demultiplexed_data/*.trimmed_renamed.fasta > demultiplexed_data/merged_seqs.fasta
../software/usearch/usearch -derep_fulllength demultiplexed_data/merged_seqs.fasta -fastaout demultiplexed_data/seqs_unique.fasta -sizeout
```

**Step5: Sort by abundance and discard singletons**
Assuming that singleton sequences (only appear once in the dataset) are likely artifacts, it is advisable to remove them from the dataset
```
../software/usearch/usearch -sortbysize demultiplexed_data/seqs_unique.fasta -fastaout demultiplexed_data/seqs_unique_sorted.fasta -minsize 2
```

**Step6: OTU clustering**
One of the most important steps in pre-processing of amplicon data is OTU clustering (or OTU 'picking'). This step consists on grouping together sequences based on similarity generating taxonomic units that roughly correspond to the species (~97% sequence identity) or genus level (~95%). During OTU clustering, the algorithm automatically *de novo* identifies chimera in the dataset. 
```
mkdir otu_result
../software/usearch/usearch -cluster_otus demultiplexed_data/seqs_unique_sorted.fasta -otus otu_result/otus.fa -uparseout otu_result/uparse.txt -relabel OTU_
```
This generates a set of representative sequences for each OTU or cluster (which are stored in the file otus.fa). Then reads (prior de-replication and either before or after quality filtering) can be mapped back onto the OTU representatives to obtain read counts for each OTU and generate an OTU table (see bellow).

We can check the size of this set:
```
grep -c '^>' otu_result/otus.fa
```


**Step7: Generate the OTU Table**
Next, all reads are mapped back onto the relatively small set of OTU representative sequences by using the command:
```
../software/usearch/usearch -usearch_global demultiplexed_data/merged_seqs.fasta -db otu_result/otus.fa -strand plus -id 0.97 -uc otu_result/read_mapping.uc
```
The output of this tool is in USEARCH cluster format (UC), a more description of UC format can be found here: https://www.drive5.com/usearch/manual/opt_uc.html
```
head otu_result/read_mapping.uc
```
and needs to be converted for downstream analyses, for example to text format, which can be loaded into R. This can be done with the following python script:
```
python ../software/usearch/uc2otutab.py otu_result/read_mapping.uc > otu_result/otu_table.txt
```
Now we can check the otu table:
```
head otu_result/otu_table.txt
```

----
### Pipeline for Reads Denoising
Reads denoising can provide highest resolution for amplicon data. This method uses exact sequence variants (100% similarity) to profile microbial composition rather than cluster reads. The most important thing is to distinguish the true reads and erroneous reads produced during PCR amplification and sequencing. Different algorithms use different statistical models to calculate the probability that one read is errorneous produced from a more abundant read. Currently, three algorithms are available for amplicon reads denoising, DADA2, Deblur and Unoise. There are relavant plugins for DADA2 and Deblur in QIIME2.

**Activate QIIME2 again**
```
conda activate qiime2-2019.7
```
Preprocess of the data from raw reads to demultiplexed reads follows the same procedure from *step1 to step2* introduced in *Pipeline for OTU clustering* 
Here we just use the demultiplexed data (data_demux-details.qza) from step2 of *Pipeline for OTU clustering* to denoise the sequencing reads. DADA2 is a prevalent algorithm to correct the sequence errors in amplicon data. .you can run it using the command:
```
mkdir asv_result
qiime dada2 denoise-paired --i-demultiplexed-seqs data_demux.qza --o-table asv_result/asv_table.qza --o-representative-sequences asv_result/asv-seqs.qza --o-denoising-stats asv_result/asv-stats.qza --p-n-threads 1 --p-trunc-len-f 0 --p-trunc-len-r 0
```
It generates three files including the asv table (similar to OTU table), representative asv sequences (similar to OTU sequences) and a file recording denoising stats. All three files are in qza format. To convert them into human-readable format, you can run commands:
```
qiime tools export --input-path asv_result/asv_table.qza --output-path asv_result/
biom convert  -i asv_result/feature-table.biom --to-tsv -o asv_result/asv_table.tsv
qiime tools export --input-path asv_result/asv-seqs.qza --output-path asv_result/
```
Now, you can check the asv table and representative asv sequences using the commands:
```
cd asv_result/
head asv_table.tsv
head dna-sequences.fasta
```
You can find the names of each asv are randomly coded, but the appearence order of each asv in both asv_table.tsv and dna-sequences.fasta is the same. So you can rename the asv with following commands:
```
awk 'BEGIN{count=1}{if(NR>2){$1="ASV_"count; print $0; count++} else{print $0}}' asv_table.tsv > asv_table_rename.tsv
awk 'BEGIN{count=1}{if(/>/){$0=">ASV_"count; print $0; count++} else{print $0}}' dna-sequences.fasta > dna-sequences_rename.fasta
```
Check the renamed files:
```
head asv_table_rename.tsv
head dna-sequences_rename.fasta
cd ..
```
----

----
# Synthetic Communities
**Pipeline for Synthetic Commuties**
The steps from raw reads to demutiplexed merged reads are the same with the step1-2 in *Pipeline for OTU clustering*. Once we get the demultiplexed samples, we can align the amplicon reads in each sample to the reference sequence of each strain with 100% identity. The abundance of each strain is calculated by counting the number of reads aligned to the reference sequnce from that strain.
Align the amplicon reads to the reference sequences with usearch:
```
mkdir syncom_result
../software/usearch/usearch -uparse_ref demultiplexed_data/merged_seqs.fasta -db reference_data/syncom_reference.fna -strand plus -uparseout syncom_result/uparse.up -threads 1
```
Modify the format of output file into uc format:
```
cat syncom_result/uparse.up |awk '{if ($2=="perfect"){print $1, $5}}' |sed 's/ /\t/g;s/^/H\t\t\t\t\t\t\t\t/g' > syncom_result/otu_clustering.uc
```
Generate OTU table with following command:
*we need to deactivate current qiime2 environment, because the default python compiler in this environment is python3, but the script to generate otu tables in Usearch is written in python2*
```
conda deactivate 
python ../software/usearch/uc2otutab.py syncom_result/otu_clustering.uc > syncom_result/otu_table.txt
```
*or we can call python2 compiler directly without deactivating qiime2 environment with the command:*
```
/netscratch/common/MPIPZ_SPP_workshop/software/miniconda2/bin/python ../software/usearch/uc2otutab.py syncom_result/otu_clustering.uc > syncom_result/otu_table.txt
```
----

# Summarizing the OTU Table
This step will print a summary of the count information on a per-sample basis and the number of unique observations/OTUs/ASVs per sample.
**Convert the tsv table to biom format**
*Biom is a file format to store information in OTU table, which is usually used in qiime1*
```
biom convert -i syncom_result/otu_table.txt -o syncom_result/otu_table.biom --table-type="OTU table" --to-hdf5
```
**Summarizing sample data**
```
biom summarize-table -i syncom_result/otu_table.biom -o syncom_result/otu_table_summary.txt
```
Now you can check the output file:
```
cat syncom_result/otu_table_summary.txt
```
**Summarizing sample data qualitatively**
```
biom summarize-table -i syncom_result/otu_table.biom --qualitative -o syncom_result/otu_table_qual_summary.txt
```
Check the output file:
```
cat syncom_result/otu_table_qual_summary.txt
```
**How to summarize the ASV table?**


 [QIIME2]: <https://qiime2.org/>
 [Usearch]: <https://www.drive5.com/usearch/>