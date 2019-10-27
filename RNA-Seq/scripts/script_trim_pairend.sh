##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

#########################trim pair-end#######################
date

FASTP_PATH=/home/yniu/Biotools

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_raw_data
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/MPIPZ_SPP_workshop/RNA-Seq_clean_data

CORENUM=16

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~
## sample names
fq=($(ls | grep _p_.*.fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p -c \
                 -h ${i}.html \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
###############################################################
