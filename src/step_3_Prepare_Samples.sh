#! /bin/bash

# set -e
# set -x

mkdir -p ${RES}/INPUT_STEP_2
mkdir -p ${RES}/TEMP

for INDIV in `awk '{print $2}' ${RES}/inds.tsv | sort | uniq`
do
	echo "	->" ${INDIV}
	awk -v I=${INDIV} -v RES=${RES} '$2==I {print RES"/OUT_PROCESS/"$1"/"$2".R1.fq.gz"}' ${RES}/inds.tsv > ${RES}/TEMP/$$.R1
	gzip -dc `tr "\n" " " < ${RES}/TEMP/$$.R1` > ${RES}/INPUT_STEP_2/${INDIV}.R1.fq
	gzip -f ${RES}/INPUT_STEP_2/${INDIV}.R1.fq
done

rm ${RES}/TEMP/$$.R1

