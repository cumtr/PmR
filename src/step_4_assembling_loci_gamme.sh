#!/bin/bash

rm -rf ${RES}/STEP_4_ASSEMBLE_LOCI
mkdir ${RES}/STEP_4_ASSEMBLE_LOCI
STACKS=${RES}/STEP_4_ASSEMBLE_LOCI

mkdir -p ${RES}/TEMP


MISMATCH_LOC_IND_MAX_START=${MISMATCH_LOC_IND_START}
MISMATCH_LOC_IND_MAX_END=${MISMATCH_LOC_IND_END}

ls -l ${RES}/STEP_3_PREPARE_SAMPLE/*.R1.fq.gz | awk '{print $9}' | awk -v seed=$RANDOM 'BEGIN{srand(seed);}{print rand()" "$0}' | sort | head -n ${NB_INDIV_M} | awk '{print $2}' >> ${RES}/TEMP/$$.ind


while ((${MISMATCH_LOC_IND_MAX_START} <= ${MISMATCH_LOC_IND_MAX_END}))
do
	mkdir -p ${STACKS}/${MISMATCH_LOC_IND_MAX_START}
	i="0"
	echo "	M = " ${MISMATCH_LOC_IND_MAX_START}
	echo Ind nb_unique_stacks nb_merged_stacks mean_cov > ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/Stat_ustacks_M${MISMATCH_LOC_IND_MAX_START}.log

	cat ${RES}/TEMP/$$.ind | while read forward
		do
        	${USTACKS} -p $THREADS -t gzfastq -f ${forward} -o ${STACKS}/${MISMATCH_LOC_IND_MAX_START} -r -i ${i} -M ${MISMATCH_LOC_IND_MAX_START} -m ${COVERAGE_LOC_MIN} 2> ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/ustacks_err
        	i=$[${i}+1]
			awk '/Assembled/' ${OUTPATH}/ustacks_err | awk '{print $2,$5}' | tr -d "\n" | awk '{print $1,$3}' > ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/nbstacks
	        	awk '/Final coverage/ {print $3}' ${OUTPATH}/ustacks_err | sed -e 's/mean=//g' | sed -e 's/;//g' > ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/cov
			chem=${forward/.AliConcat.fq/}
			name=${chem/${RES}/STEP_3_PREPARE_SAMPLE\//}
			echo ${name} `cat ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/nbstacks` `cat ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/cov` >> ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/Stat_ustacks_M${MISMATCH_LOC_IND_MAX_START}.log
		done


	rm ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/cov
	rm ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/nbstacks
	rm ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/ustacks_err
	MISMATCH_LOC_IND_MAX_START=$[${MISMATCH_LOC_IND_MAX_START}+1]
done

rm ${RES}/TEMP/$$.ind


MISMATCH_LOC_IND_MAX_START=${MISMATCH_LOC_IND_START}
MISMATCH_LOC_IND_MAX_END=${MISMATCH_LOC_IND_END}

if [ -f ${STACKS}/Mismatch_allowed.txt ] ; then
    rm ${STACKS}/Mismatch_allowed.txt
fi

touch ${STACKS}/Mismatch_allowed.txt

while ((${MISMATCH_LOC_IND_MAX_START} <= ${MISMATCH_LOC_IND_MAX_END}))
        do
        awk '{print $3}' ${STACKS}/${MISMATCH_LOC_IND_MAX_START}/Stat_ustacks_M${MISMATCH_LOC_IND_MAX_START}.log > ${STACKS}/$$
	paste -d "\t" ${STACKS}/Mismatch_allowed.txt ${STACKS}/$$ > $$.1
	mv $$.1 ${STACKS}/Mismatch_allowed.txt
	rm ${STACKS}/$$

	MISMATCH_LOC_IND_MAX_START=$[${MISMATCH_LOC_IND_MAX_START}+1]

	done


R --vanilla <<EOF

        TAB = read.table(paste0("${STACKS}","/Mismatch_allowed.txt"), skip = 1)
        TAB_Part = apply(TAB, 1, function(x){1-(x/x[1])})

        pdf(paste0("${STACKS}","/M_choice.OutStep4.pdf") , width = 10, height = 5)

        boxplot(t(TAB_Part), border = "grey50", names = c(${MISMATCH_LOC_IND_START}:${MISMATCH_LOC_IND_END}), main = "Part of polymorphics stacks\naccording to M values (R1)")
        lines(apply(TAB_Part, 1, function(x){median(x)}), col = "blue", lwd = 2)

        dev.off()
EOF

