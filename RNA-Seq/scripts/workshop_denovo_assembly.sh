##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

date

WORKSHOP_PATH=/netscratch/common/MPIPZ_SPP_workshop

TRINITY_PATH=/opt/share/software/bin
SELECT_PATH=${WORKSHOP_PATH}/RNA-Seq/RNA-Seq_select_data

ASSEMBLY_PATH=~/trinity_assembly_select

if [ ! -d "${SELECT_PATH}" ]; then
    mkdir ${SELECT_PATH}
fi

CORENUM=2

cd ${SELECT_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~de novo assembly~~~~~~~~~~~~~~~~~~~~~~~~~
${TRINITY_PATH}/Trinity --seqType fq \
               --left Col0_p_rep1_select_R1.fq.gz,flg22_p_rep1_select_R1.fq.gz \
               --right Col0_p_rep1_select_R2.fq.gz,flg22_p_rep1_select_R2.fq.gz \
               --output ${ASSEMBLY_PATH} \
               --trimmomatic \
               --normalize_by_read_set \
               --min_contig_length 400 \
               --max_memory 20G \
               --CPU ${CORENUM}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
