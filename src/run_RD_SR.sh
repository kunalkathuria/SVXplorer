#!/bin/bash
WORK_DIR=${32}

#Uncomment this to extract split-read file if not available, with appropriate path replacement
#samtools view -h ${14} \
#    | ~/store/lumpy-sv-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
#    | samtools view -Sb - \
#    > /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#
#samtools sort -n -@32 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.ns.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#rm /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam

cp $WORK_DIR/ClassifiedVariantMap.txt $WORK_DIR/VariantMap.txt
cp $WORK_DIR/All_Variants_O.txt $WORK_DIR/All_Variants.txt
python add_SR_hash.2.py $1 $2 $3 $4 ${13} ${12} ${10} ${18} ${21} $WORK_DIR
cp $WORK_DIR/VariantMap_SR.txt $WORK_DIR/VariantMap.txt
cp $WORK_DIR/All_Variants_SR.txt $WORK_DIR/All_Variants.txt
python SetCover_mq.py $6 ${11} ${15} ${16} ${31} $WORK_DIR
python WriteBed.o2.py $WORK_DIR/DisjSetCover_S.txt $WORK_DIR

if [ -d $WORK_DIR/sr_results ]
then
        rm -r $WORK_DIR/sr_results
fi
mkdir $WORK_DIR/sr_results
mv $WORK_DIR/*.bedpe $WORK_DIR/sr_results

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
	done < <(tr -d '\r' < $WORK_DIR/bam_stats.txt)
	RD_SLOP=$((meanIL+DISC_D ))
	echo $RD_SLOP
	cat $WORK_DIR/All_Variants.txt | grep -v INV | grep -v known > $WORK_DIR/inDels.txt # all INDELS
	sort -k6,6 $WORK_DIR/inDels.txt > $WORK_DIR/inDels_S.txt
	#./formatRD.sh #make first 2 vars in here user-input
	python add_RD.o2.py $5 $8 ${15} ${16} ${17} ${19} $RD_SLOP ${24} ${23}
	cp $WORK_DIR/VariantMap_RD.txt $WORK_DIR/VariantMap.txt
	cp $WORK_DIR/All_Variants_RD.txt $WORK_DIR/All_Variants.txt
fi

if [ $9 -eq 1 ]
then
	python SetCover_mq.py $6 ${11} ${15} ${16} ${31} $WORK_DIR
else
	python DisjointSetCover.py $WORK_DIR
fi

python WriteBed.o8.o.py $WORK_DIR/DisjSetCover_S.txt ${14} ${22} ${23} ${25} ${26} ${27} ${28} ${29} ${30} $WORK_DIR
echo "WB Done"
cat $WORK_DIR/deletions.bedpe $WORK_DIR/tandemDuplications.bedpe $WORK_DIR/inversions.bedpe $WORK_DIR/insertions.bedpe $WORK_DIR/unknowns.bedpe > $WORK_DIR/All_SVs
python ../scripts/bed2vcf.py $WORK_DIR/All_SVs
cp $WORK_DIR/VariantMap_O.txt $WORK_DIR/VariantMap.txt
cp $WORK_DIR/All_Variants_O.txt $WORK_DIR/All_Variants.txt

