##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##


#########################download SRA################################
date

BIN_PATH=/usr/bin
SRATOOLS_PATH=/extDisk1/Biotools/sratoolkit.2.10.0-centos_linux64/bin

RAW_PATH=/extDisk1/RESEARCH/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_raw

CORENUM=10

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~single-end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ref: https://www.nature.com/articles/nature21417

SRARaw=("SRR4295742" "SRR4295743" "SRR4295744"
        "SRR4295748" "SRR4295749" "SRR4295750"
        "SRR4295754" "SRR4295755" "SRR4295756"
        "SRR4295778" "SRR4295779" "SRR4295790"
        "SRR4295784" "SRR4295785" "SRR4295786"
        "SRR4295790" "SRR4295791" "SRR4295792")

SRAName=("Col0_s_rep1_batch1" "Col0_s_rep2_batch1" "Col0_s_rep3_batch1"
         "MeJA_s_rep1_batch1" "MeJA_s_rep2_batch1" "MeJA_s_rep3_batch1"
         "flg22_s_rep1_batch1" "flg22_s_rep2_batch1" "flg22_s_rep3_batch1"
         "Col0_s_rep1_batch2" "Col0_s_rep2_batch2" "Col0_s_rep3_batch2"
         "MeJA_s_rep1_batch2" "MeJA_s_rep2_batch2" "MeJA_s_rep3_batch2"
         "flg22_s_rep1_batch2" "flg22_s_rep2_batch2" "flg22_s_rep3_batch2")

for i in {1..9}; do

    ${SRATOOLS_PATH}/fasterq-dump ${SRARaw[i]} -e 12 \
                    -O ${RAW_PATH}

    ${BIN_PATH}/gzip ${SRARaw[i]}.fastq

    ${BIN_PATH}/mv ${SRARaw[i]}.fastq.gz ${SRAName[i]}.fq.gz

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~trim~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date

FASTP_PATH=/home/yniu/Biotools

RAW_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_raw
CLEAN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_clean

CORENUM=16

cd ${RAW_PATH}

## sample names
fq=($(ls | grep _s_.*.fq.gz))

for i in ${fq[@]}; do

    echo "Trimming ${i}."
    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p \
                 -h ${i%%.*}.html \
                 -i ${i} \
                 -o ${CLEAN_PATH}/${i}

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~alignment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/ath
CLEAN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_clean
ALIGN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_batch/RNA-Seq_batch_align

KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.45.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
RM_PATH=/bin/rm
MOVE_PATH=/bin/mv

CORENUM=45
SPECIES='ath'

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~single end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep _s_.*.fq.gz))
fqnames=($(echo "${fq[@]%%.*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    --single -l 200 -s 20\
                    ${CLEAN_PATH}/${i}.fq.gz
done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
