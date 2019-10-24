##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

date

FASTP_PATH=/netscratch/common/MPIPZ_SPP_workshop/software/fastp

RAW_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_smallraw_data
CLEAN_PATH=~/RNA-Seq_clean_data

if [ ! -d "${CLEAN_PATH}" ]; then
    mkdir ${CLEAN_PATH}
fi

CORENUM=2

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~single-end~~~~~~~~~~~~~~~~~~
## sample names
fq=($(ls | grep _s_.*.fq.gz))

for i in ${fq[@]}; do

    echo "Trimming ${i}."
    ${FASTP_PATH}/fastp -w ${CORENUM} \
                     -z 6 \
                     -p \
                     -h ${CLEAN_PATH}/${i%%.*}.html \
                     -i ${i} \
                     -o ${CLEAN_PATH}/${i}

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~pair-end~~~~~~~~~~~~~~~~~~
## sample names
fq=($(ls | grep _p_.*.fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."
    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p -c \
                 -h ${CLEAN_PATH}/${i}.html \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
