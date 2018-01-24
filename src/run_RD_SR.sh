#!/usr/bin/env bash
set -eu
WORK_DIR=${28}
SRC=${29}
VERBOSE=${31}
THREADS=${32}
#Uncomment this to extract split-read file if not available, with appropriate path replacement
#samtools view -h ${12} \
#    | ~/store/lumpy-sv-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
#    | samtools view -Sb - \
#    > /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#
#samtools sort -n -@ $THREADS -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.ns.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
#rm /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/splitters.us.bam
echo $WORK_DIR, $SRC, $VERBOSE
$SRC/addSplitReads.py -s $2 -f $3 -m ${11} -q ${10} -c ${8} -i ${16} -t ${19} $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/variantMap.pe.txt $WORK_DIR/SVC_debug/allVariants.pe.txt $1
python $SRC/uniqueSuppFilter.py -g $9 -k $5 $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/allDiscordants.txt
#intermediate output here

#if [ -d $WORK_DIR/SVC_results/sr_results ]
#then
#        rm -r $WORK_DIR/SVC_results/sr_results
#fi
#mkdir $WORK_DIR/SVC_results/sr_results
#mv $WORK_DIR/SVC_results/*.bedpe $WORK_DIR/SVC_results/sr_results

if [ $6 -eq 1 ]
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
	python $SRC/addRD.py $4 $7 ${13} ${14} ${15} $RD_SLOP ${18} ${22} ${21}
fi
python $SRC/uniqueSuppFilter.py -g ${9} -k $5 $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/allDiscordants.txt

python $SRC/covPUFilter.py -a ${20} -b ${21} -c ${26} -e ${25} -f ${30} -g ${27} -v $VERBOSE \
    $WORK_DIR/SVC_results $WORK_DIR/SVC_debug/variants.uniqueFilter.txt \
    $WORK_DIR/SVC_debug/allVariants.pe_sr.txt $WORK_DIR/SVC_debug/variantMap.pe_sr.txt ${12} $WORK_DIR/SVC_debug/bamStats.txt
echo "WB Done"

cat $WORK_DIR/SVC_results/deletions.bedpe $WORK_DIR/SVC_results/tandemDuplications.bedpe $WORK_DIR/SVC_results/inversions.bedpe $WORK_DIR/SVC_results/insertions.bedpe $WORK_DIR/SVC_results/unknowns.bedpe > $WORK_DIR/SVC_results/All_SVs
python $SRC/../scripts/bed2vcf.py $WORK_DIR/SVC_results/All_SVs
