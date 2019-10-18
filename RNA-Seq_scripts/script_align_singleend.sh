date

REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/ath
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_clean_data
ALIGN_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_align_data

KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.45.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
RM_PATH=/bin/rm
MOVE_PATH=/bin/mv

CORENUM=45
SPECIES='ath'

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athht2index/genome \
                  -U ${CLEAN_PATH}/${i}.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
