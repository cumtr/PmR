STEP#!/bin/bash

mkdir -p ${RES}/STEP_2_DEMULTIPLEX

for forward in ${RES}/STEP_1_FITLER_READS/*/*R1.cut.fastq.gz
do
	reverse=${forward/R1.cut.fastq.gz/R2.cut.fastq.gz}
	RUN=${forward/_R1.cut.fastq.gz/}
	LIB=`basename ${RUN}`
        echo ""
        echo ""
	echo "------------- Demultiplexing " ${LIB} " -------------"
        echo ""
        echo ""
 	mkdir -p ${RES}/STEP_2_DEMULTIPLEX/${LIB}
 	${PROCESS_RADTAG} -1 ${forward} -2 ${reverse} -D -o ${RES}/STEP_2_DEMULTIPLEX/${LIB} -b ${DATA}/${LIB}.barcode --renz_1 ${ENZYME_1} --renz_2 ${ENZYME_2} -c -r -i 'gzfastq' -t ${MAX_LENGTH}
	#mv ${RES}/STEP_2_DEMULTIPLEX/${LIB}/process_radtags.log ${RES}/STEP_2_DEMULTIPLEX/${LIB}/process_radtags.${LIB}.log
	awk 'BEGIN{f=FILENAME"."nb++;}/^$/{f=FILENAME"."nb++;next;}{print > f;}' ${RES}/STEP_2_DEMULTIPLEX/${LIB}/process_radtags.${LIB}.log

R --vanilla <<EOF
	
	TAB = read.table("${RES}/STEP_2_DEMULTIPLEX/${LIB}/process_radtags.${LIB}.log.3", skip = 1)
	TAB_Order = TAB[order(TAB[,6]),]
	TAB_Reduced = rbind(TAB_Order[,6] , TAB_Order[,3]-TAB_Order[,6])

	pdf("${RES}/STEP_2_DEMULTIPLEX/${LIB}/${LIB}.log.3.pdf", width = 10, height = 7)
	layout(matrix(c(1,1,1,1,1,1,1,1,1,2,3,4),ncol = 4))
	x = barplot(height = TAB_Reduced, col = c("blue", "white"), main = paste0("Number of reads for each individual\n","${LIB}"))
	text(cex=0.8, x=x+0.5, y=0-4, TAB_Order[,2], xpd=TRUE, srt=45, pos=2)
	legend(x = 2, y = max(TAB_Order[,3]), legend = c("Total", "Retained"), fill = c("white","blue"), bty = "n")
	
	plot.new()
	pie(rowSums(TAB_Reduced), col = c("blue","white"), radius = 1, init.angle = 90, labels = round(c(rowSums(TAB_Reduced)/sum(rowSums(TAB_Reduced)))*100, digits = 2), main = "Mean prop.")
	dev.off()

EOF
	echo "	Done"

done

##############################################################################

##############################################################################

echo ""
echo ""
echo "	Concatenate R1 and Rem1"

for forward in `find ${RES}/STEP_2_DEMULTIPLEX -type f -name \*.1.fq.gz`
do
    echo $forward
    name=${forward/1.fq.gz/}
    if [[ ! $forward =~ rem ]]
        then
        cat $name*1* > ${name}R1.fq.gz
    fi
done

echo "  Done"
echo ""
echo ""

##############################################################################

##############################################################################


for barcode in ${DATA}/*.barcode
   do
        LIB=`basename ${barcode}`
        LIB=${LIB/.barcode/}
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${RES}/${LIB}.OutStep2.pdf ${RES}/STEP_2_DEMULTIPLEX/${LIB}/${LIB}.log.3.pdf
   done


find ${RES}/STEP_2_DEMULTIPLEX -type f -name \*.R1.fq.gz | awk -F'/' 'BEGIN{OFS="\t"}{id=$NF;gsub(/.R1.fq.gz$/, "", id); print $(NF-1), id, id;"$NF;"}' | sort > ${RES}/inds.tsv



