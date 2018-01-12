#!/bin/bash


OUTPATH=${RES}/STEP_5_7_RUN_STACKS/
mkdir -p ${RES}/TEMP

for forward in ${OUTPATH}/*.snps.tsv.gz
    do
        name=${forward/.snps.tsv.gz}
        echo -s ${name} >> ${RES}/TEMP/list
    done

${CSTACKS} -b 1 -o ${OUTPATH} -p ${THREADS} -n ${MISMATCH_CATALOG_MAX} `cat ${RES}/TEMP/list`

rm ${RES}/TEMP/list
