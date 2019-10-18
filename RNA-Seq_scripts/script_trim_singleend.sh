#########################trim single-end#######################
date

FASTP_PATH=/home/yniu/Biotools

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_raw_data
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_clean_data

CORENUM=16

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~
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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
###############################################################
