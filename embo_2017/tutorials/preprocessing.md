
### Hands-on session - Pre-processing of raw amplicon data

#### Formatting of barcodes
Depending on the instrument, barcode sequences might be provided as a separate
FASTQ file or directly in the (forward and / or reverse) read headers. It is
possible to convert between these formats, e.g. by using the QIIME command:

``add_barcode_to_label.py ./data/001_forward_reads.fastq.gz ./data/001_barcodes.fastq.gz ./forward_reads.fastq``

``add_barcode_to_label.py ./data/001_reverse_reads.fastq.gz ./data/001_barcodes.fastq.gz ./reverse_reads.fastq``

Now the barcodes given in the file 001_barcodes.fastq.gz, e.g.

``zcat ../data/001_barcodes.fastq.gz | head``

have been appended to the read headers

``head forward_reads.fastq``

Converting between these formats can be useful as some scripts or toolkits
will only take one or the other as input.

#### Merge (assemble) paired-end reads

The first step in pre-processing of Illumina paired-end data consists in
joining forward and reverse reads. In most (but not all) cases, there will
be an overlap which will allow to correct errors in the terminal (low quality)
positions of the reads. This can be problematic with genetic loci with
intervals of variable length between forward and reverse primers (e.g. ITS)
but it is generally not an issue for bacterial 16S genes.

It is possible to merge paired-end reads using the following command:

``usearch -fastq_mergepairs forward_reads.fastq -reverse reverse_reads.fastq -fastqout joined.fastq``

After joining, the barcodes can be again extracted from the read headers.
This is only necessary when using different toolkits to avoid compatibility
issues

``extract_barcodes.py -f joined.fastq -c barcode_in_label --char_delineator 'BC=' --bc1_len 12 -o ./``

#### Sample demultiplexing

Using the information that matches sample IDs with their corresponding barcodes
(file ../data/001_mapping.txt), it is possible to assign each read to a sample
(and discard those erroneous barcodes), e.g. by using the following QIIME script:

``split_libraries_fastq.py -i joined.fastq -b barcodes.fastq --rev_comp_mapping_barcodes --max_barcode_errors 2 -o ./ -m ./data/001_mapping.txt -q 30 --phred_offset 33``

Note that while it is possible to generate a separate FASTQ or FASTA file for
each sample, which can be useful to submit the raw data to a public repository
such as ENA, here we append all sequences to the same file (joined.fastq) for
convenience, as these will all be jointly processed in subsequent steps.

For compatibility between QIIME and USEARCH, we need to relabel the sequence
headers, and can use this simple script:

``./scripts/relabel_ids.sh seqs.fna seqs.fasta``

To see some statistics on the demultiplexing (and in this case also quality
filtering) step, one can check the log file

``cat split_library_log.txt``

#### Dereplication

Sequences without errors that come from the same biological entity are identical.
Ideally, most of the dataset would consist of many replicates of the same few
(real) template sequences. Although this is unfortunately not the case (~50%
of the reads only, due to errors and artifacts), it is possible to substantially
speed-up downstream processing by removing replicates (while keeping count),
e.g. by using the following command in USEARCH:

``usearch -derep_fulllength seqs.fasta -fastaout seqs_unique.fasta -sizeout``

We can see the reduction in size of these FASTA files by running:

``grep -c ">" seqs.fasta``

``grep -c ">" seqs_unique.fasta``

#### Sort by abundance and discard singletons

Assuming that singleton sequences (only appear once in the dataset) are likely
artifacts, it is advisable to remove them from the dataset, e.g. by using the
following command:

``usearch -sortbysize seqs_unique.fasta -fastaout seqs_unique_sorted.fasta -minsize 2``

and likewise check the reduction in file size

``grep -c ">" seqs_unique_sorted.fasta``

#### OTU clustering

One of the most important steps in pre-processing of amplicon data is OTU
clustering (or OTU 'picking'). This step consists on grouping together sequences
based on similarity generating taxonomic units that roughly correspond to
the species (~97% sequence identity) or genus level (~95%).

We can perform OTU clustering using UPARSE by running the command:

``usearch -cluster_otus seqs_unique_sorted.fasta -otus otus_uparse.fasta -id 0.97``

This generates a set of representative sequences for each OTU or cluster (which
are stored in the file otus_uparse.fasta). After some denoising and filtering
of this representative sequences, reads (prior de-replication and either before
or after quality filtering) can be mapped back onto the OTU representatives to
obtain read counts for each OTU and generate an OTU table (see bellow).

We can check the size of this set, which at this stage is likely to include
many artifacts, e.g. from chimeric sequences:

``grep -c ">" otus_uparse.fasta``

#### OTU de-noising and chimera detection

Chimeras (artifactual sequences formed from two or more biological sequences
joined together during PCR) are a common problem with amplicon sequencing data.

Detecting chimeric sequences can be computationally demanding, so we generally
avoid doing this on the complete set of input sequences. Instead, it is enough
to filter chimeric OTUs from the dataset, which will automatically remove all
chimeric reads from the final OTU table.

There are multiple methods for chimera detection, e.g.: ChimeraSlayer
(Haas *et al.*, 2011), DECIPHER (Wright *et al.*, 2012), Perseus
(Quince *et al.*, 2011) or UCHIME (Edgar *et al.*, 2011). As an example,
we can apply UCHIME by using the following command:

``usearch -uchime_ref otus_uparse.fasta -db ./data/cs_gold.fa -strand plus -nonchimeras otus_nc_uparse.fasta -threads 1``

and check the number of OTU representative sequences flagged as chimeric and
removed from the dataset

``grep -c ">" otus_nc_uparse.fasta``

There are more steps that can further remove artifacts and errors from the
set of OTUs, for example by aligning representatives to public databases and
removing all sequences that do not match anything previously sequenced above
a certain threshold of sequence identity (we generally use 75% but less
stringent values might also be meaningful).

We next rename the final set of representative OTU sequences, e.g. directly
in bash using 'awk':

``cat otus_nc_uparse.fasta | awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' >> rep_seqs.fasta``

``grep ">" rep_seqs.fasta``

#### Generate the (count) OTU table

Next, all reads are mapped back onto the relatively small set of OTU representative
sequences by using the command

``usearch -usearch_global seqs.fasta -db rep_seqs.fasta -strand plus -id 0.97 -uc read_mapping.uc``

The output of this tool is in USEARCH cluster format (UC)

``head read_mapping.uc``

and needs to be converted for downstream analyses, for example to text format,
which can be loaded into R. This can be done with the following python script:

``python ./scripts/usearch/uc2otutab.py read_mapping.uc > otu_table.txt``

Some toolkits (e.g. QIIME) require OTU tables in biom format, which are more
efficient than plain text to encode large matrices of count data. These files
can be generated by typing:

``biom convert -i otu_table.txt -o otu_table.biom --table-type="OTU table" --to-json``

