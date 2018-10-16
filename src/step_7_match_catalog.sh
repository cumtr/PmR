#!/bin/bash

OUTPATH=${RES}STEP_5_7_RUN_STACKS/

for forward in ${OUTPATH}/*.snps.tsv.gz
do
    if [[ ! $forward =~ catalog ]]
    then
        forw=${forward/.snps.tsv.gz/}
        ${SSTACKS} -p ${THREADS} -c ${OUTPATH} -s ${forw} -o ${OUTPATH}/
    fi
done


tsv2bam -P ${OUTPATH} -M ${POP_INFOS} -t ${THREADS} 

gstacks -P ${OUTPATH} -M ${POP_INFOS} -t ${THREADS} 
