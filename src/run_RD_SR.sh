#!/bin/bash
WORK_DIR=${30}
SRC=${31}
VERBOSE=${33}
echo $WORK_DIR, $SRC, $VERBOSE
#Uncomment this to extract split-read file if not available, with appropriate path replacement
#samtools view -h ${14} \
#    | ~/store/lumpy-sv-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
#    | samtools view -Sb - \
#    > /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#
#samtools sort -n -@32 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.ns.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#rm /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam

$SRC/addSplitReads.py -r $2 -s $3 -f $4 -m ${13} -q ${12} -c ${10} -i ${18} -t ${21} $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/variantMap.pe.txt $WORK_DIR/SVC_debug/allVariants.pe.txt $1
python $SRC/uniqueSuppFilter.py -g ${11} -k $6 $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/allDiscordants.txt
#intermediate output here

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
	python $SRC/add_RD.o2.py $5 $8 ${15} ${16} ${17} ${19} $RD_SLOP ${24} ${23}
fi
python $SRC/uniqueSuppFilter.py -g ${11} -k $6 $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/allDiscordants.txt
python $SRC/covPUFilter.py -a ${22} -b ${23} -c ${28} -e ${27} -f ${32} -g ${29} -v $VERBOSE \
    $WORK_DIR/SVC_results $WORK_DIR/SVC_debug/variants.uniqueFilter.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt ${14} $WORK_DIR/SVC_debug/bamStats.txt
echo "WB Done"
cat $WORK_DIR/SVC_results/deletions.bedpe $WORK_DIR/SVC_results/tandemDuplications.bedpe $WORK_DIR/SVC_results/inversions.bedpe $WORK_DIR/SVC_results/insertions.bedpe $WORK_DIR/SVC_results/unknowns.bedpe > $WORK_DIR/SVC_results/All_SVs
python $SRC/../scripts/bed2vcf.py $WORK_DIR/SVC_results/All_SVs
