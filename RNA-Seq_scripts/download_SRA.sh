#########################download SRA################################
date

BIN_PATH=/usr/bin
SRATOOLS_PATH=/extDisk1/Biotools/sratoolkit.2.10.0-centos_linux64/bin

RAW_PATH=/extDisk1/RESEARCH/MPIPZ_SPP_workshop/RNA_Seq_raw_data/

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~single-end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ref: https://www.nature.com/articles/nature21417

SRARaw=("SRR4295742" "SRR4295743" "SRR4295744"
        "SRR4295748" "SRR4295749" "SRR4295750"
        "SRR4295754" "SRR4295755" "SRR4295756")

SRAName=("Col0_s_rep1" "Col0_s_rep2" "Col0_s_rep3"
         "MeJA_s_rep1" "MeJA_s_rep2" "MeJA_s_rep3"
         "flg22_s_rep1" "flg22_s_rep2" "flg22_s_rep3")

for i in {0..8}; do

    ${SRATOOLS_PATH}/fasterq-dump ${SRARaw[i]} -e 12 \
                    -O ${RAW_PATH}

    ${BIN_PATH}/gzip ${SRARaw[i]}.fastq

    ${BIN_PATH}/mv ${SRARaw[i]}.fastq.gz ${SRAName[i]}.fq.gz

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~pair-end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ref: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1007499

SRARaw=("SRR7865799" "SRR7865800" "SRR7865801"
        "SRR7865805" "SRR7865806" "SRR7865807")

SRAName=("Col0_p_rep1" "Col0_p_rep2" "Col0_p_rep3"
         "flg22_p_rep1" "flg22_p_rep2" "flg22_p_rep3")

for i in {0..5}; do
    ${SRATOOLS_PATH}/fasterq-dump ${SRARaw[i]} -e 12 \
                    -O ${RAW_PATH}

    # ${BIN_PATH}/gzip ${SRARaw[i]}.fastq

    # ${BIN_PATH}/mv ${SRARaw[i]}.fastq.gz ${SRAName[i]}.fq.gz
done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
#####################################################################
