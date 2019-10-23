date

WORKSHOP_PATH=/netscratch/common/MPIPZ_SPP_workshop

KALLISTO_PATH=${WORKSHOP_PATH}/software/kallisto-0.46.0
HISAT2_PATH=${WORKSHOP_PATH}/software/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
BIN_PATH=/bin

REF_PATH=${WORKSHOP_PATH}/RNA-Seq/RNA-Seq_index
CLEAN_PATH=~/RNA-Seq_clean_data
ALIGN_PATH=~/RNA-Seq_align_data

if [ ! -d "${ALIGN_PATH}" ]; then
    mkdir ${ALIGN_PATH}
fi

CORENUM=2
SPECIES='ath'

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pair end~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep _p_.*.fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    ${CLEAN_PATH}/${i}_R1.fq.gz ${CLEAN_PATH}/${i}_R2.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athht2index/genome \
                  -1 ${CLEAN_PATH}/${i}_R1.fq.gz \
                  -2 ${CLEAN_PATH}/${i}_R2.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${BIN_PATH}/rm ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
