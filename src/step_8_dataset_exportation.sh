#!/bin/bash

INPUT=${RES}/STEP_5_7_RUN_STACKS/
OUTPUT=${RES}/STEP_8_FINAL_OUTPUT/

mkdir -p ${OUTPUT}

STACKS=${RES}/STACKS
INP_REP=${STACKS}/CATALOG_${M_CHOOSEN}/


if [ -e  ${POP_INFOS} ]
then
	${POPULATIONS} -b 1 -t ${THREADS} -O ${OUTPUT}  -P ${INPUT} -M ${POP_INFOS} -p ${NB_POP} -r ${PROP_POP} -a ${MAF} --write_random_snp --vcf --structure --genepop --phylip --fstats --plink --phylip_var_all
else
	${POPULATIONS} -b 1 -t ${THREADS} -O ${OUTPUT}  -P ${INPUT} -r ${PROP_POP} -a ${MAF} --write_random_snp --vcf --structure --genepop --phylip --fstats --plink --phylip_var_all
fi

sed -n '/Distribution of population loci\./,/Distribution of confounded loci at catalog locus\./p' ${OUTPUT}/batch_1.populations.log | awk '! /#/' > ${OUTPUT}/batch_1.populations.txt
sed -n '/Distribution of population loci after applying locus constraints\./,/Distribution of confounded loci at catalog locus\./p' ${OUTPUT}/batch_1.populations.log | awk '! /#/' > ${OUTPUT}/batch_1.populations.retained.txt

R --vanilla <<EOF

  table_nbloci<-read.table("${OUTPUT}/batch_1.populations.txt", header=F)
  table<-data.frame(1,1)
  for (i in 0:(nrow(table_nbloci)-1))
	{table[i+1,1]<-table_nbloci[i+1,1]
	 table[i+1,2]<-sum(table_nbloci[(nrow(table_nbloci)-i):nrow(table_nbloci),2])}
  table_nbloci2<-read.table("${OUTPUT}/batch_1.populations.retained.txt", header=F)
  table2<-table
  table2[,2]<-0
  for (i in 0:(nrow(table_nbloci2)-1))
        {table2[i+1,2]<-sum(table_nbloci2[(nrow(table_nbloci2)-i):nrow(table_nbloci2),2])}
  pdf(paste0("${OUTPUT}","R1_Rarefaction_curve.pdf"), width=(5+round(nrow(table_nbloci2)/5)))
    barplot(rev(table[,2]), col="blue", names.arg=table[,1],cex.names=0.5)
    barplot(rev(table2[,2]), col="white", add=TRUE)
  dev.off()

EOF
