 ## EMBO practical course - Plant microbiota
 ### March / April 2017 | MPIPZ Cologne

This repository contains the materials used for the computational and theoretical
sessions taught at the 1st EMBO practical course on Plant Microbiota in March /
April 2017 at the Max Planck Institute for Plant Breeding Research in Cologne,
Germany.

This is a private repository made available to the participants of the course
for their personal use and training. Please, contact
[me](mailto:garridoo@mpipz.mpg.de) or any or the other course organizers in
case you wish to share these materials with other people. You should feel free
to reuse or adapt any or all of the code found here without any guarantees from
the author(s). If you do so in a publication, a mention to the author(s) in the
acknowledgements section would be greatly appreciated. Lastly a note of caution:
if you are not entirely familiar with the basic concepts of the analysis that
you intend to perform it is generally a very good idea to revisit the reference
material or to approach a qualified computational biologist for advice; these
tutorials are meant to be a general guide and NOT a step-by-step detailed
protocol to use without careful thought.

If you find that some of the materials are missing, or if you find any mistakes
or inconsistencies, please, open up an 'issue' at the repository instead of
sending [me](mailto:garridoo@mpipz.mpg.de) an email, unless urgent. You can
do this on the github repository website by clicking on the 'issue' icon on
the top bar. Then, I will try to fix the problem as soon as possible.
Note also that I will be adding still missing presentations and materials
from external trainers over the next few days / weeks.

Bellow is a brief outline of the contents of this repository with direct
links to the markdown tutorials, PDF or PPT slides and presentations, scripts
and example data used during the duration of the course.

#### Contents of this repository

* [tutorials](https://github.com/garridoo/embo_2017/tree/master/tutorials):
contains the markdown files of the step-by-step tutorials.
* [slides](https://github.com/garridoo/embo_2017/tree/master/slides):
contains the PDF or PPT presentations of the theoretical sessions.
* [data](https://github.com/garridoo/embo_2017/tree/master/data):
in this folder you should download and extract the example data (see bellow).
* [scripts](https://github.com/garridoo/embo_2017/tree/master/scripts):
contains the scripts necessary for processing, plotting and statistical tests.
* [figures](https://github.com/garridoo/embo_2017/tree/master/figures):
empty folder which plotting scripts use as an output for figures.

#### Wetlab documentation

The slides used for the lectures related to the wetlab part of the course
can be found in a single PDF file
[here](https://github.com/garridoo/embo_2017/blob/master/slides/wetlab.pdf).
In the same folder you can find the protocols of the
[first](https://github.com/garridoo/embo_2017/blob/master/slides/protocols_1.pdf)
and
[second weeks](https://github.com/garridoo/embo_2017/blob/master/slides/wetlab.pdf).

#### Downloading the example data

First, you should download a tarball containing all the example data (too large
to host in github) from our server. The data is hosted
[here](http://www.at-sphere.com/downloads/embo_data_2017.tar.gz).
If you are using a linux terminal in a remote server, you can type:

`wget http://www.at-sphere.com/downloads/embo_data_2017.tar.gz`

you can extract the contents of this file into your data folder by using the
following command:

`tar -zxvf embo_data_2017.tar.gz`

Note that this toy data was generated from both published and unpublished
data and doesn't necessarily has any biological meaning. It is meant to be
used only to test the scripts provided bellow. I recommend you directly try
with your own dataset.

#### Introduction to amplicon data processing

The slides of the introductory lecture "[State-of-the-art approaches for
amplicon data processing](https://github.com/garridoo/embo_2017/blob/master/slides/pre_processing_slides.pdf)"
can be found [here](https://github.com/garridoo/embo_2017/blob/master/slides/pre_processing_slides.pdf).
These slides contain a general overview of the most relevant steps required
to analyze marker gene data, including raw data, pre-processing, sequence
analyses and diversity / statistical tests.

#### Raw data pre-processing

The first tutorial correspond to the steps required to pre-process raw
amplicon data and generation of the OTU tables. This includes sequence
merging (joining pair-end reads), quality filtering, demultiplexing,
OTU clustering, de-noising, chimera removal, etc. The step-by-step
markdown tutorial can be found
[here](https://github.com/garridoo/embo_2017/blob/master/tutorials/pre_processing.md).

#### Taxonomic classification of marker gene data

Once we obtain a OTU table and the corresponding representative sequences,
we can assign to each OTU a taxonomic label, using a variety of approaches.
The slides from the lecture
"[Taxonomic characterization and diversity analyses](https://github.com/garridoo/embo_2017/blob/master/slides/taxonomy_slides.pdf)"
contain a brief summary of the most widely used methods. The markdown tutorial
with example commands and explanations can be found
[here](https://github.com/garridoo/embo_2017/blob/master/tutorials/taxonomy.md).

#### Alpha-and beta diversity analyses

One of the key analyses that can be performed on community data are the
so-called diversity analyses. The basics concepts are summarized in the talk
"[Alpha- and beta-diversity](https://github.com/garridoo/embo_2017/blob/master/slides/diversity.pdf)".
Conducting these analysis begins by calculating the corresponding alpha
(within sample) and beta (between sample) diversity estimates. Using
OTU tables (generally normalized or rarefied) we can obtain these indices
e.g. by using the steps described in the corresponding tutorial, which can be
found
[here](https://github.com/garridoo/embo_2017/blob/master/tutorials/diversity.md).

#### Visualizing community data

Using the OTU tables and the diversity indices calculated before, we can
visualize community data, e.g. by plotting
[relative abundances using stacked barplots](https://github.com/garridoo/embo_2017/blob/master/scripts/stacked_barplots.R),
[alpha-diversity boxplots](https://github.com/garridoo/embo_2017/blob/master/scripts/alpha_diversity.R) or
[PCoA / MNDS scatter plots of beta-diversity](https://github.com/garridoo/embo_2017/blob/master/scripts/beta_diversity.R), etc.
An overview of the basic concepts behind multidimensional scaling and ordination
methods can be found
[here](https://github.com/garridoo/embo_2017/blob/master/slides/ordination_methods.pdf).
Other examples of other visualization tools such as ternary graphs, Manhattan
plots, phylogenies, etc. can be found in the
[code repositories of previous publications](http://www.mpipz.mpg.de/R_scripts)
of the department.

#### Whole-genome assembly and annotation

An alternative way to explore the diversity of a microbial community is to
analyze the whole genomic sequences of isolated community members. The slides
of the talk
"[Overview of bioinformatics tools for whole-genome sequence analysis](https://github.com/garridoo/embo_2017/blob/master/slides/wgs.pdf)"
contain a brief overview of the basic steps involved in bacterial genome
assembly and annotation, and phylogenomics. The corresponding tutorial to
this section can be found
[here](https://github.com/garridoo/embo_2017/blob/master/tutorials/wgs.md).

#### Multidimensional scaling and ordination methods

An overview of the basic concepts behind multidimensional scaling and ordination
methods can be found
[here](https://github.com/garridoo/embo_2017/blob/master/slides/ordination_methods.pdf).
Some examples of how to assess the contribution of environmental variables to
diversity (e.g. using constrained ordination approaches) can be found
[here](https://github.com/garridoo/embo_2017/blob/master/scripts/cpcoa.R).

#### Reproducibility and replicability

Reproducibility in computational biology is a topic of heated debate.
To ensure maximum transparency, it is necessary to adhere to good practices
which are at present unfortunately not required for publication in most
journals. These include depositing code and scripts in public repositories
that could be used by other researchers to replicate your results.
[Here](https://github.com/garridoo/embo_2017/blob/master/slides/reproducibility.pdf)
are some slides where we briefly discussed this important topic.
