##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

BIN_PATH=/usr/bin
SAMTOOLS_PATH=/home/yniu/Biotools/samtools_1.9/bin
SEQTK_PATH=/home/yniu/Biotools/seqtk

CLEAN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_clean_data
ALIGN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_align_data
SELECT_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_select_data

## FLS2 -- 5:18791575-18795908
## MKK3 -- 5:16181628-16184827
## MEKK1 -- 4:5403671-5407320

if [ ! -d "${SELECT_PATH}" ]; then
    mkdir ${SELECT_PATH}
fi

cd ${SELECT_PATH}

## extract reads names from bam
${SAMTOOLS_PATH}/samtools view -b ${ALIGN_PATH}/Col0_p_rep1_ath_hisat2.bam 5:18791575-18795908 4:5403671-5407320 > ${SELECT_PATH}/Col0_p_rep1_ath_hisat2_select.bam

${SAMTOOLS_PATH}/samtools view -b ${ALIGN_PATH}/flg22_p_rep1_ath_hisat2.bam 5:18791575-18795908 4:5403671-5407320 > ${SELECT_PATH}/flg22_p_rep1_ath_hisat2_select.bam

## extract reads from fastq
${SAMTOOLS_PATH}/samtools view Col0_p_rep1_ath_hisat2_select.bam \
    | ${BIN_PATH}/awk '{print $1}' \
    | ${BIN_PATH}/sort \
    | ${BIN_PATH}/uniq > Col0_p_rep1_rnames.txt

${SAMTOOLS_PATH}/samtools view flg22_p_rep1_ath_hisat2_select.bam \
    | ${BIN_PATH}/awk '{print $1}' \
    | ${BIN_PATH}/sort \
    | ${BIN_PATH}/uniq > flg22_p_rep1_rnames.txt

${SEQTK_PATH}/seqtk subseq ${CLEAN_PATH}/Col0_p_rep1_R1.fq.gz Col0_p_rep1_rnames.txt \
    | ${BIN_PATH}/awk '{print (NR%4 == 1) ? "@61DFRAAXX100204:1:100:10494:" ++i "/1" : $0}' \
    | /bin/gzip > ${SELECT_PATH}/Col0_p_rep1_select_R1.fq.gz

${SEQTK_PATH}/seqtk subseq ${CLEAN_PATH}/Col0_p_rep1_R2.fq.gz Col0_p_rep1_rnames.txt \
    | ${BIN_PATH}/awk '{print (NR%4 == 1) ? "@61DFRAAXX100204:1:100:10494:" ++i "/2" : $0}' \
    | /bin/gzip > ${SELECT_PATH}/Col0_p_rep1_select_R2.fq.gz

${SEQTK_PATH}/seqtk subseq ${CLEAN_PATH}/flg22_p_rep1_R1.fq.gz flg22_p_rep1_rnames.txt \
    | ${BIN_PATH}/awk '{print (NR%4 == 1) ? "@61DFRAAXX100204:1:100:10495:" ++i "/1" : $0}' \
    | /bin/gzip > ${SELECT_PATH}/flg22_p_rep1_select_R1.fq.gz

${SEQTK_PATH}/seqtk subseq ${CLEAN_PATH}/flg22_p_rep1_R2.fq.gz flg22_p_rep1_rnames.txt \
    | ${BIN_PATH}/awk '{print (NR%4 == 1) ? "@61DFRAAXX100204:1:100:10495:" ++i "/2" : $0}' \
    | /bin/gzip > ${SELECT_PATH}/flg22_p_rep1_select_R2.fq.gz

date

