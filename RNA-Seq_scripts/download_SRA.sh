#########################download SRA################################
date

BIN_PATH=/usr/bin
SRATOOLS_PATH=/extDisk1/Biotools/sratoolkit.2.10.0-centos_linux64/bin

RAW_PATH=/extDisk1/RESEARCH/MPIPZ_SPP_workshop/RNA_Seq_raw_data/

cd ${RAW_PATH}

SARRaw=("SRR4295742" "SRR4295743" "SRR4295744"
        "SRR4295748" "SRR4295749" "SRR4295750"
        "SRR4295754" "SRR4295755" "SRR4295756")

SARName=("Col0_rep1" "Col0_rep2" "Col0_rep3"
         "MeJA_rep1" "MeJA_rep2" "MeJA_rep3"
         "flg22_rep1" "flg22_rep2" "flg22_rep3")

for i in {0..8}; do

    ${SRATOOLS_PATH}/fasterq-dump ${SARRaw[i]} -e 12 \
                    -O ${RAW_PATH}

    ${BIN_PATH}/gzip ${SARRaw[i]}.fastq

    ${BIN_PATH}/mv ${SARRaw[i]}.fastq.gz ${SARName[i]}.fq.gz

done

date
#####################################################################
