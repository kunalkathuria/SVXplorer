#!/bin/bash
WORK_DIR=${32}
SRC=${33}

#Uncomment this to extract split-read file if not available, with appropriate path replacement
#samtools view -h ${14} \
#    | ~/store/lumpy-sv-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
#    | samtools view -Sb - \
#    > /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#
#samtools sort -n -@32 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.ns.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#rm /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam

python $SRC/add_SR_hash.2.py $1 $2 $3 $4 ${13} ${12} ${10} ${18} ${21} $WORK_DIR/SVC_debug
python $SRC/SetCover_mq.py $6 ${11} ${15} ${16} ${31} $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/variantMap.pe_sr.txt $WORK_DIR/SVC_debug/allVariants.pe_sr.txt
python $SRC/WriteBed.o2.py $WORK_DIR/SVC_debug/variants.uniqueFilter.txt $WORK_DIR/SVC_results $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/allVariants.pe_sr.txt

if [ -d $WORK_DIR/SVC_results/sr_results ]
then
        rm -r $WORK_DIR/SVC_results/sr_results
fi
mkdir $WORK_DIR/SVC_results/sr_results
mv $WORK_DIR/SVC_results/*.bedpe $WORK_DIR/SVC_results/sr_results

if [ $7 -eq 1 ]
then
	counter=0
	while read -r line;
	do
        	counter=$((counter+1))
        	if [[ $counter -eq 2 ]]; then
                	echo $line
                	meanIL=${line%.*}
        	elif [[ $counter -eq 6 ]]; then
                	echo ${line%.*}
                	DISC_D=${line%.*}
        	fi
	done < <(tr -d '\r' < $WORK_DIR/SVC_debug/bamStats.txt)
	RD_SLOP=$((meanIL+DISC_D ))
	echo $RD_SLOP
	cat $WORK_DIR/SVC_debug/allVariants.pe_sr.txt | grep -v INV | grep -v known > $WORK_DIR/SVC_debug/inDels.txt # all INDELS
	sort -k6,6 $WORK_DIR/SVC_debug/inDels.txt > $WORK_DIR/SVC_debug/inDels_S.txt
	#./formatRD.sh #make first 2 vars in here user-input
	python $SRC/add_RD.o2.py $5 $8 ${15} ${16} ${17} ${19} $RD_SLOP ${24} ${23}
fi

if [ $9 -eq 1 ]
then
	python $SRC/SetCover_mq.py $6 ${11} ${15} ${16} ${31} $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/variantMap.pe_sr.txt $WORK_DIR/SVC_debug/allVariants.pe_sr.txt
else
	python $SRC/DisjointSetCover.py $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/variantMap.pe_sr.txt $WORK_DIR/SVC_debug/allVariants.pe_sr.txt
fi

python $SRC/WriteBed.o8.o.py $WORK_DIR/SVC_debug/variants.uniqueFilter.txt ${14} ${22} ${23} ${25} ${26} ${27} ${28} ${29} ${30} $WORK_DIR/SVC_results $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt ${34}
echo "WB Done"
cat $WORK_DIR/SVC_results/deletions.bedpe $WORK_DIR/SVC_results/tandemDuplications.bedpe $WORK_DIR/SVC_results/inversions.bedpe $WORK_DIR/SVC_results/insertions.bedpe $WORK_DIR/SVC_results/unknowns.bedpe > $WORK_DIR/SVC_results/All_SVs
python $SRC/../scripts/bed2vcf.py $WORK_DIR/SVC_results/All_SVs

