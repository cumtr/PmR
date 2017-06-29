#!/bin/bash

#set -e
#set -x

OUTPATH=${RES}/STACKS_OUTPUT/

mkdir -p ${OUTPATH}

MISMATCH=${M_CHOOSEN}

rm -f ${OUTPATH}/Stat_ustacks_M${MISMATCH}.log
rm -f ${OUTPATH}/ustacks_err
rm -f ${OUTPATH}/nbstacks
rm -f ${OUTPATH}/cov

i="0"
echo Ind nb_unique_stacks nb_merged_stacks mean_cov > ${OUTPATH}/Stat_ustacks_M${MISMATCH}.log

for forward in ${RES}/INPUT_STEP_2/*.R1.fq.gz
do
	echo `basename ${forward}`
        ${USTACKS} -p $THREADS -t gzfastq -f ${forward} -o ${OUTPATH} -d -i ${i} -M ${MISMATCH} -m ${COVERAGE_LOC_MIN} 2> ${OUTPATH}/ustacks_err
        i=$[${i}+1]
        awk '/stacks merged into/ { print $1,$5}' ${OUTPATH}/ustacks_err  > ${OUTPATH}/nbstacks
        awk '/After merging/ { print $6}' ${OUTPATH}/ustacks_err | sed -e 's/;//g'  > ${OUTPATH}/cov
        chem=${forward/.R1.fq.gz/}
        name=${chem/${OUT_PROCESS}\//}
        echo ${name} `cat ${OUTPATH}/nbstacks` `cat ${OUTPATH}/cov` >> ${OUTPATH}/Stat_ustacks_M${MISMATCH}.log
done

R --vanilla <<EOF

	TABLE = read.table("${OUTPATH}/Stat_ustacks_M${MISMATCH}.log", header = T)
	TABLE = TABLE[order(TABLE[,3]),]
	names=unlist(lapply(strsplit(as.character(TABLE[,1]), "INPUT_STEP_2/", fixed=TRUE), function(x) x[2]))

 	pdf(paste0("${RES}","Indiv_FragCover.OutStep5.pdf"), width=(5+round(nrow(TABLE)/5)), height=10)

	par(mar=c(8,6,4,6))
	b = barplot(TABLE[,2], col = "white", names.arg = names, las = 3, cex.names = 0.7, axes = F)
	axis(2,las = 2)
	barplot(TABLE[,3], col = "blue", add = T, axes = F)
	mtext("Nb Stacks", side=2, line=3, cex.lab=1,las=3)
	
	par(new=TRUE)
	plot(1,type = "n", xlim = c(min(b), max(b)+1), ylim = c(0, max(TABLE[,4])), axes = F, xlab ="", ylab="")
	c = b+0.5
	points(TABLE[,4]~(c), pch = 20, type = "b")
	axis(4,las = 2)
	mtext("Median cov", side=4, line=3, cex.lab=1, las=3)
	dev.off()

EOF


rm ${OUTPATH}/cov
rm ${OUTPATH}/nbstacks
rm ${OUTPATH}/ustacks_err


