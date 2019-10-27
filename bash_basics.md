# Introduction to Shell

## What is a Shell ?

A shell is a command-line interpreter which enables the user to directly access the file structure and execute commands on a computer or server without using a Grafical User Interface (GUI)

All commands are entered in a terminal in textformat.


## 1 Basic Commands

We start with some simple commands that are used to do basic operations in the bash. Like moving around or show contents of files and folders.

Commands are usually structured the following way:

**command** \[ *(-positional) a* \] \[ *--optional b* \] \[-(-) flag\]

```shell
ls -l -h some_folder
```

With certainty you will have to add at least one required (sometimes called positional) argument. Further you might have the option to give optional arguments (e.g. choose another algorithm) or set flags which can turn certain options on and off.

So in our example above our command is **ls**. It lists information about files. Our required argument is the folder we want to list. Here: *some\_folder*. And we set two flags **-l** will format the files as a list and **-h** will output the file sizes in a human readable format.

### Finding help

To find out what exactly the options for a command are or what additional possibilities you have you can usually try one of the following things:

- ###### Help

  try your command with -h, --h, -help or --help

- ###### Linux manpages

  man \[command\]

- ###### Google


### Exercise 1.1

Open your terminal and try the ls command with increasing complexity.

```Shell
# ls in the current directory
ls

# ls in the root folder and list
ls / -l

# ls in the root folder list and file sizes human readable
ls / -l -h
```

#### Basic Commands Overview

###### Moving Around and Creating Files and  Directories
- **ls** *list directory*
- **cd** *command directory*
- **touch** *make a new empty file*
- **mkdir** *make a new empty directory*
- **rm** *remove file / folder*

###### Show File Content
- **cat** *print file to the standard output*
- **cut** *extract columns from a file*
- **head** *print only the first x lines of a file*
- **tail** *print only the tail of a file*
- **less** *interactive file browser*
- **wc** *count words, lines*

###### Connect to a Server
- **ssh** *connect to a server with your username and password*

###### Working with compressed Files

- **gzip** *compress files*
- **gunzip** *decompress files*
- **zcat** *like cat but for compressed files*
- **zless** *like less but for compressed files*

#### Good to know:

A "**.**" denotes the current directory and "**..**" devotes the directory which is the parent of the current one in the file structure. Further "**~**" is your home directoy and "**/**" is the root directoy of the system.

E.g.: You want to copy a file to the current directory.

```shell
# cp <source> <destination>

# copy distant file to current directoy
cp some_other_folder/some_other_file .

# copy file from current directoy to parent folder
cp some_file ..
```
Another important symbol is the asterisk "**\***" which can be seen as a wildcard when selecting files.

Lets say you want to select all file of a directory with the ending ".faa" then you simply write "\*.faa"

So if you would like to delete all these files you would simply write:

```Shell
# remove all files with the ending .faa
rm *.faa
```

### Let's connect to the server

As a shell can be an interface for any computer - even distant ones - we now want to log in to a server of the institute.

In general the ssh command works like this:

ssh username@servername

```Shell
#connect to a server called dell-node-13 in the institute spp01 only example login name

ssh ssp01@dell-node-13

# next you will be asked for your password. Type it in and press enter.

```

You are now on the server! If you look around you will notice, that there is no data in your home directory. In the next step we will copy some data in our home directory.


### Exercise 1.2

The data for the workshop is all stored in the common folder */netscratch/common/MPIPZ_SPP_workshop/*. We will now copy the data folder for the shell workshop into your home directory.

What does the **-r** flag indicate? Use ```man cp``` to find out.

```shell
# copy data folder
cp -r /netscratch/common/MPIPZ_SPP_workshop/shell_data shell_workshop

# goto data directory
cd shell_workshop
```

### Exercise 1.3

Make a new directory with the name *temporary*. Touch 2 new files (named bacteria.faa and bacteria.gff) in the folder. Delete the bacteria.faa file first and then the complete temporary folder.

```shell
#Make the directory
mkdir temporary

# Touch new files
touch temporary/bacteria.fasta
touch temporary/bacteria.gff

# File created by piping output of one command into a new file
echo -P ">gene1\nACTTATAGGGA\n" > temporary/another_bacteria.fasta

# remove files and directories
rm temporary/bacteria.gff
rm -r temporary
```

### Exercise 1.4

Cat cut head and tail are all commands to stream a file to the standard output ( screen ). Where cat will put out the whole content of a file, head will give you the first n lines and tail the n last lines. If you have a file which is made out of columns you can use cut to choose specific columns to put out.

```Shell
# show contens of files
cat fasta_files/genome_1.fasta

# show head or tail of files
head -n 50 fasta_files/genome_1.fasta
tail -n 5 fasta_files/genome_1.fasta

# choose specific columns
cut -f 2 blast_results/Blast10_24.txt
 ```

### Exercise 1.5
Sometimes files can be rather large and you might want to quickly browse a file. So instead of opening it in a GUI text editor you could use less for example.

Try open a file with less. To quit less you can press **q**.

You can also use less to search in a file. Simply type **/** followed by your search term.

```shell
#open a file with less
less fasta_files/genome_1.fasta

# press q to quit
# press / + search term to search
```

### Exercise 1.6

To quickly check the number of words or lines in a file you can use *wc*.

```shell
# count number of words
wc -w fasta_files/genome_1.fasta

#count number of lines
wc -l fasta_files/genome_1.fasta

```

### Exercise 1.7

Besides the standard shell commands there are a lot of other programs that you can use.

For you it might be interesting to start an interactive R session or run an R script

```shell
# run R
R
# run an R script
Rscript scripts/example.R

# alternativly run R source("scripts/example.R")


# run python and python script
python

python scripts/hello_world.py

# and now with your name as an argument
python scripts/hello_world.py your_name
```

### Exercise 1.8

Use gzip und gunzip to compress and decompress a fasta file

``` shell
# compress
gzip fasta_files/genome_1.fasta

# decompress
gunzip fasta_files/genome_1.fasta.gz

# compress but keep original file
gzip -c fasta_files/genome_1.fasta > fasta_files/genome_1_compressed.gz

# compare sizes with:
 ls -lh fasta_files

# try cat with compressed file
cat fasta_files/genome_1_compressed.gz

# now
zcat fasta_files/genome_1_compressed.gz
```

# 2 Advanced File Manipulations

In this chapter you will learn advanced file manipulations.

#### Good to know:

On the shell you can concatenate the output of several commands by piping the output of one command into the other. For that purpose one uses piping symbols:

- the pipe symbol **|**  is used to direct the output of one command to another
- **>** this symbol is used to direct the output of a command to a text file. But be careful all existing contents of a file will be deleted.
- **>>** is similar to **>** the difference is that it will be appended to the file it is directed to.

These piping symbols are the real power of the shell. Because you can combine many commands and achieve operations very quick on many files.

#### Advanced File Manipulators Command Overview

###### File Stream Manipulators
- **grep** *use regular expressions to filter certain lines*
- **awk** *filter and manipulate lines*
- **sed** *manipulate lines*
- **tr** *trim lines, replace and delete characters*
- **sort** *sort lines of a file stream*
- **uniq** *merge uniqe elements of a stream*

###### Text Editors
- **vim** *text editor which is versatile in use, but a bit harder to learn*
- **nano** *text editor which is easier to use but limited in applications*

### Exercise 2.1

The **grep** command is a nice tool to search files for strings or search terms. You can search by string or a regular expressions (RE). The latter is a special language to formulate patterns that you want to look for. For example if you want to look for all lines that start with a number the RE would be: "^\d" or "^[0-9]". Where "^" stands for beginning of the line and "\d" or "[0-9]" for any number. For more information please have a look in the appendix.

In the first exercise you are going to get all sequence identifiers from a fasta file. As you know sequence identifiers in a fasta file always start with a ">". All you have to do is grep for this symbol.

```Shell
# grep the headers
grep ">" fasta_files/genome_1.fasta

# count sequences in fasta file
grep -c ">" fasta_files/genome_1.fasta

# alternative:
grep ">" fasta_files/genome_1.fasta  | wc -l

# write identifiers to files
grep ">" fasta_files/genome_1.fasta  > documents/identifier_list.txt
```

Another example for using grep is looking for motifs. You can simply grep a fasta file motiv. In the following example we will look for the TATA box motif.

```Shell
# we start by greping the motif tata and look at the first 50 lines of the output
# the option i stands for case insensitive
# the color option is to see wich subsequences match your search term
grep -i "TATA" --color=always fasta_files/genome_1.fasta | head -n 50

# this gives us already a nice overview but usually the complete tata box motif will be in a form like TATA(A ot T)A(A or T). Grep is well suited for such fuzzy search terms.
# Here we use further option n and o. Where n will put out the line number and o will put out only the sequence that matched your search term.
# If you use these square brackets in a regular expresions they mean any of the symbol it contains. So in the below case the sequence "TATA" must be followed by an ["A" or an "T"]

grep -ino --color=always "TATA[AT]A[AT]" fasta_files/genome_1.fasta | head -n 50

#In the previous example count how many motifs you have found
#Write 2000 lines of the output from the second example into a text file called tata_with_line.txt. Move it to the documents folder.
#Do the same without line numbers to a file called tata.txt
#Try to modify the above expression such that any nucleotide can follow the "TATA" motif.
```

### Exercise 2.2

Another nice tool to manipulate files is **awk**. It reads the lines of an input file and spits them into fields. You can then work with these fields and apply filter operations, arithmetic operations and conditional expressions.

**awk** has a special syntax that you have to use. As said above each line gets split into fields. If you want to refer to the whole line you use **$0** and the fields are then **$1 ... $n** where n is the number of fields

For our purposes the following command structure is important to learn

```Shell
# basic awk
awk '{ what to do }' input_file

# sometimes you need to change the field separator for this u can use the -F options. In this example it is explicitly the tab separator '\t'

awk -F '\t' '{ what to do }' input_file
```

Now lets do a very simple example. We only want to select certain columns from our file. And in the second step we will only print columns that fullfill a certain condition

```shell
# print field 1, 2 and 4
awk '{ print $1"\t"$2"\t"$4 }' blast_results/Blast10_24.txt

# print field 1, 2 and 3 if the value in the third column is bigger than 90
awk '{ if ( $3 > 90) {print $1"\t"$2"\t"$3} }' blast_results/Blast10_24.txt

# print field 1,2,3 and 12(bit score) where the bitscore is bigger then 100. Save the result to blast_results/results_bitscore_gt100.txt
```

### Exercise 2.3

The last file stream manipulator we will learn is **sed**. It can be used for a variety of file stream opertations but the most important feature is to find and replace strings in a line.

```Shell
# basic sed command will replace the first occurence of search string
sed 's/search pattern/replace with/'

# for replacing all occurences of the string you can use g for global
sed 's/search pattern/replace with/g'

# if you want to make the changes inplace in the file you can use the -i (in place) option but you should only do this if you are sure about what you are doing.
sed -i 's/search pattern/replace with/g'

```

Now try it with some of the example files.

```Shell
# replace all spaces with tabs

sed 's/\t/:::/g' blast_results/results_bitscore_gt100.txt

# replace

sed 's/$/\tnew colum/g' blast_results/results_bitscore_gt100.txt

```

### Exercise 2.4
**sort** and **uniq** can be quite helpful in many situations. Sort is for sorting text files and uniq collapses unique lines

```shell
# sort the filtered blast results by percentage identity
# ascending
sort -k 3 blast_results/results_bitscore_gt100.txt
# descending
sort -k 3 -r  blast_results/results_bitscore_gt100.txt

# count the number of uniq tata box motifs in the tata.txt file. For this we have to sort the file first and then use unique with the count option
sort documents/tata.txt | uniq -c  
```
## Mounting a server directory to your laptop

```shell
# install sshfs - exit ssh first!
sudo apt install sshfs

# make directory in your homefolder(!) [cd ~]
mkdir server_home

# mount server home to folder server_home
sshfs spp01@dell-node-13 server_home

# to unmount the directory later on do:
sudo umount server_home

```



## Regular Expressions

##### Some important symbols:
- **+** expression is there at least once or infinitely many times
- **\*** expression not there or infinitely many times
- **{n}** expression is there n times
- **[ ]** any of the symbols in the square brackets
- **^** start of the line
- **$** end of the line
