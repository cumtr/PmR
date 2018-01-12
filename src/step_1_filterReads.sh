#!/bin/bash


mkdir -p ${RES}/STEP_1_FITLER_READS/

for f in ${DATA}/*_R1.fastq.gz 
do
    lib=${f/_R1.fastq.gz/}
    outname=`basename $lib`
    echo "	processing " $lib " ....."
    
    mkdir -p ${RES}/STEP_1_FITLER_READS/${outname}
	
    ${BBDUK} in1=${lib}_R1.fastq.gz in2=${lib}_R2.fastq.gz out1=${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R1.bbduk.fastq.gz out2=${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R2.bbduk.fastq.gz ref=${ADAP} ktrim=r k=20 mink=10 hdist=1 tpe tbo overwrite=true bhist=${RES}/STEP_1_FITLER_READS/${outname}/bhist_${outname}.txt qhist=${RES}/STEP_1_FITLER_READS/${outname}/qhist_${outname}.txt aqhist=${RES}/STEP_1_FITLER_READS/${outname}/aqhist_${outname}.txt lhist=${RES}/STEP_1_FITLER_READS/${outname}/lhist_${outname}.txt >& ${RES}/STEP_1_FITLER_READS/${outname}/${outname}.bbduk.log
    
    mv ${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R1.bbduk.fastq.gz ${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R1.cut.fastq.gz
    mv ${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R2.bbduk.fastq.gz ${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R2.cut.fastq.gz
    echo "	Done" 

echo ${RES}/STEP_1_FITLER_READS/${outname}

gzip -dc ${RES}/STEP_1_FITLER_READS/${outname}/${outname}_R1.cut.fastq.gz | awk '(FNR%4)==2{l[length($0)]+=1}END{for (i in l) print i, l[i];}' > ${RES}/STEP_1_FITLER_READS/${outname}/reads.length.txt

LENGTH=${RES}/STEP_1_FITLER_READS/${outname}/reads.length.txt


R --vanilla <<EOF

  table_compo<-read.table("$LENGTH", header=F, row.names=1)
  tutu<-as.character(row.names(table_compo)[order(as.numeric(row.names(table_compo)))])
  table_compo <- as.data.frame(table_compo[order(as.numeric(row.names(table_compo))),])
  row.names(table_compo)<-tutu  

  pdf("${RES}/STEP_1_FITLER_READS/${outname}/Reads_length.pdf", width = 10, height = 7)

  par(mar=c(8, 6, 4, 6))
 
  b = barplot(table_compo[,1], col = "blue", xlab="Reads length after adapters trimming", names.arg = row.names(table_compo), las = 3, cex.names=0.5, axes = F)
  axis(2, las=2)
  mtext("Nb reads", side = 2, line = 4, cex.lab = 1, las = 3)
 
  toto<-NULL
  for (i in as.numeric(row.names(table_compo)))
  {toto<-c(toto, sum(table_compo[as.numeric(row.names(table_compo))<i,1]))}
  tutu<-1-(toto/sum(table_compo[,1]))

  par(new=TRUE)
  plot(1, type="n", xlim = c(min(b), max(b)), ylim = c(0, 1), axes = F, xlab="", ylab="")
  c = b + 0.5
  points(tutu~(c), pch = 20, type = "p", cex = 0.7, col = "grey70")
  points(tutu~(c), pch = 20, type = "l", cex = 0.7, col = "grey70")
  axis(4, las = 2)
  mtext("Percent of retained reads", side=4, line=3, cex.lab=1, las=3)

  abline(v=b[which(as.numeric(row.names(table_compo))==$MAX_LENGTH)], col = "grey70", lty = 2)
  text((b[which(as.numeric(row.names(table_compo))==$MAX_LENGTH)]+1.5),1 ,paste0(as.character(round(tutu[which(as.numeric(row.names(table_compo))==$MAX_LENGTH)], digit = 2)*100)," %"), col = "grey70")

  dev.off()

EOF


done

##########################################################

##########################################################

echo "Processing quality statistics on reads"


for forward in ${RES}/STEP_1_FITLER_READS/*/*_R1.cut.fastq.gz
   do
	lib=${forward/_R1.cut.fastq.gz/}
	lib=`basename ${lib}`
	COMPO=${RES}/STEP_1_FITLER_READS/${lib}/bhist_${lib}.txt
 	QUAL=${RES}/STEP_1_FITLER_READS/${lib}/qhist_${lib}.txt
 

R --vanilla <<EOF


  table_compo<-read.table("$COMPO", header=F, row.names=1)
  prop<-t(as.matrix(table_compo[c(1:(nrow(table_compo)/2)),]))
  table_quality<-read.table("$QUAL", header=F, row.names=1)
  pdf(paste0("${RES}/STEP_1_FITLER_READS/${lib}/","${lib}","_quality_nucleotide.pdf"), width=10, height=15)
  layout(matrix(c(1:2), ncol = 1))
  plot(table_quality[,1]~as.numeric(row.names(table_quality)),xlim=c(0,(nrow(table_quality)+3)),ylim=c(0,(max(table_quality[,1])+10)),col="blue", pch=20, xlab="Nucleotide position along the read", ylab="Quality score", main="${lib}")
  ## Nucleotid proportion at each reads position
  barplot(prop, col=c("blue", "white", "green", "pink", "black"), xlab="nucleotide composition along the read")
  legend("topright", fill=c("blue", "white", "green", "pink", "black"), legend=c("A","C","G","T","N"))
  dev.off()


EOF

done


for barcode in ${DATA}/*.barcode
   do
        LIB=`basename ${barcode}`
        LIB=${LIB/.barcode/}
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${RES}/${LIB}.OutStep1.pdf ${RES}/STEP_1_FITLER_READS/${LIB}/Reads_length.pdf ${RES}/STEP_1_FITLER_READS/${LIB}/${LIB}_quality_nucleotide.pdf 
   done
