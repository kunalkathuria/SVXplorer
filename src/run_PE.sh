#!/bin/bash
set -eux

MINCS=$1
WORK_DIR=${13}
SRC=${14}

echo "Classification Start"
cat $WORK_DIR/SVC_debug/allClusters.txt | awk '$2 >= '$MINCS'' > \
    $WORK_DIR/SVC_debug/allClusters.thresh.txt
sort -k4,4 -k5,5n $WORK_DIR/SVC_debug/allClusters.thresh.txt > \
    $WORK_DIR/SVC_debug/allClusters.ls.txt
sort -k7,7 -k8,8n $WORK_DIR/SVC_debug/allClusters.thresh.txt > \
    $WORK_DIR/SVC_debug/allClusters.rs.txt
time ($SRC/consolidatePEClusters.py -r $2 -s $3 -f $4 -v ${15} $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/allClusters.ls.txt \
    $WORK_DIR/SVC_debug/allClusters.rs.txt $WORK_DIR/SVC_debug/clusterMap.txt)

python $SRC/uniqueSuppFilter.py -g $8 -k $5 $WORK_DIR/SVC_debug \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/variantMap.pe.txt \
    $WORK_DIR/SVC_debug/allVariants.pe.txt $WORK_DIR/SVC_debug/allDiscordants.txt

python $SRC/WriteBed.o2.py $WORK_DIR/SVC_debug/variants.uniqueFilter.txt $WORK_DIR/SVC_results \
    $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/allVariants.pe.txt

if [ -d $WORK_DIR/SVC_results/pe_results ]
then
	rm -r $WORK_DIR/SVC_results/pe_results
fi
mkdir $WORK_DIR/SVC_results/pe_results
mv $WORK_DIR/SVC_results/*.bedpe $WORK_DIR/SVC_results/pe_results
