#!/bin/bash


OUTPATH=${RES}/STEP_5_7_RUN_STACKS/
mkdir -p ${RES}/TEMP

if [ -e  ${RES}/TEMP/list ]
then 
	rm ${RES}/TEMP/list
fi

for forward in ${OUTPATH}/*.snps.tsv.gz
    do
        name=${forward/.snps.tsv.gz}
        echo -s ${name} >> ${RES}/TEMP/list
    done

${CSTACKS} -o ${OUTPATH} -p ${THREADS} -n ${MISMATCH_CATALOG_MAX} `cat ${RES}/TEMP/list`

rm ${RES}/TEMP/list
