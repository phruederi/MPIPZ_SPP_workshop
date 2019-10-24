##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

date

WORKSHOP_PATH=/netscratch/common/MPIPZ_SPP_workshop

KALLISTO_PATH=${WORKSHOP_PATH}/software/kallisto-0.46.0
HISAT2_PATH=${WORKSHOP_PATH}/software/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
BIN_PATH=/bin

REF_PATH=${WORKSHOP_PATH}/RNA-Seq/RNA-Seq_index
CLEAN_PATH=${WORKSHOP_PATH}/RNA-Seq/RNA-Seq_smallclean_data
ALIGN_PATH=~/RNA-Seq_align_data

if [ ! -d "${ALIGN_PATH}" ]; then
    mkdir ${ALIGN_PATH}
fi

CORENUM=2
SPECIES='ath'

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~single end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Q1.pattern for single-end##
fq=($(ls ${CLEAN_PATH} | grep _s_.*.fq.gz))
fqnames=($(echo "${fq[@]%%.*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."

    ## Q2: command for kallisto alignment? ##
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    --single -l 200 -s 20\
                    ${CLEAN_PATH}/${i}.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athht2index/genome \
                  -U ${CLEAN_PATH}/${i}.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ## Q3: command for samtools sort reads? ##
    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${BIN_PATH}/rm ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
