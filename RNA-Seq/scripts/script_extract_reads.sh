date

SAMTOOLS_PATH=/home/yniu/Biotools/samtools_1.9/bin
ALIGN_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_align_data
SELECT_PATH=/netscratch/common/MPIPZ_SPP_workshop/RNA-Seq/RNA-Seq_select_data

## FLS2 -- 5:18791575-18795908
## MKK3 -- 5:16181628-16184827
## MEKK1 -- 4:5403671-5407320

if [ ! -d "${SELECT_PATH}" ]; then
    mkdir ${SELECT_PATH}
fi

${SAMTOOLS_PATH}/samtools view ${ALIGN_PATH}/Col0_p_rep1_ath_hisat2.bam 5:18791575-18795908 4:5403671-5407320 > ${SELECT_PATH}/Col0_p_rep1_ath_hisat2_select.bam

${SAMTOOLS_PATH}/samtools view ${ALIGN_PATH}/flg22_s_rep1_ath_hisat2.bam 5:18791575-18795908 4:5403671-5407320 > ${SELECT_PATH}/flg22_s_rep1_ath_hisat2_select.bam

date

