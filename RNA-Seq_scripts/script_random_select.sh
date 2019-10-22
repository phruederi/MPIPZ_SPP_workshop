date

BIN_PATH=/bin
SEQTK_PATH=/home/yniu/Biotools/seqtk

RAW_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_raw_data
SMALLRAW_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_smallraw_data

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~pair-end~~~~~~~~~~~~~~~~~~~~~
fq=($(ls | grep _p_.*.fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "Sampling ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${SEQTK_PATH}/seqtk sample -s12345 ${i}_R1.fq.gz 1000 > ${SMALLRAW_PATH}/${i}_small_R1.fq
    ${BIN_PATH}/gzip ${SMALLRAW_PATH}/${i}_small_R1.fq

    ${SEQTK_PATH}/seqtk sample -s12345 ${i}_R2.fq.gz 1000 > ${SMALLRAW_PATH}/${i}_small_R2.fq
    ${BIN_PATH}/gzip ${SMALLRAW_PATH}/${i}_small_R2.fq

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~single-end~~~~~~~~~~~~~~~~~~~~~
fq=($(ls | grep _s_.*.fq.gz))
fqnames=($(echo "${fq[@]%%.*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "Sampling ${i}.fq.gz"

    ${SEQTK_PATH}/seqtk sample -s12345 ${i}.fq.gz 1000 > ${SMALLRAW_PATH}/${i}_small.fq
    ${BIN_PATH}/gzip ${SMALLRAW_PATH}/${i}_small.fq

done

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
