#! /bin/bash

# set -e
# set -x

rm -rf ${RES}/STEP_3_PREPARE_SAMPLE 
mkdir ${RES}/STEP_3_PREPARE_SAMPLE
mkdir -p ${RES}/TEMP

for INDIV in `awk '{print $3}' ${RES}/inds.tsv | sort | uniq`
do
	echo "	->" ${INDIV}
	awk -v I=${INDIV} -v RES=${RES} '$3==I {print RES"/STEP_2_DEMULTIPLEX/"$1"/"$2".R1.fq.gz"}' ${RES}/inds.tsv > ${RES}/TEMP/$$.R1
	cat `tr "\n" " " < ${RES}/TEMP/$$.R1` > ${RES}/STEP_3_PREPARE_SAMPLE/${INDIV}.R1.fq.gz
	#gzip -dc `tr "\n" " " < ${RES}/TEMP/$$.R1` > ${RES}/STEP_3_PREPARE_SAMPLE/${INDIV}.R1.fq
	#gzip -f ${RES}/STEP_3_PREPARE_SAMPLE/${INDIV}.R1.fq
done

rm ${RES}/TEMP/$$.R1

