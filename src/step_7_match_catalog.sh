#!/bin/bash

OUTPATH=${RES}/STEP_5_7_RUN_STACKS/

for forward in ${OUTPATH}/*.snps.tsv.gz
do
    if [[ ! $forward =~ catalog ]]
    then
        forw=${forward/.snps.tsv.gz/}
        ${SSTACKS} -b 1 -p ${THREADS} -c ${OUTPATH}/batch_1 -s ${forw} -o ${OUTPATH}/
    fi
done

