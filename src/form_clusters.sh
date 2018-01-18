#!/bin/bash
set -eux
WORK_DIR=${16}
SRC=${17}

time ($SRC/writeDiscordantFragments.py -p $2 -a $3 -t $5 -r $4 -n $6 -i $9 -c ${10} \
    -m ${11} -v ${18} $WORK_DIR/SVC_debug $WORK_DIR/SVC_debug/aln1s.bam \
    $WORK_DIR/SVC_debug/aln2s.bam $1)
echo Sorting
sort -T $WORK_DIR/SVC_debug -k 2,2 -k 4,4 -k 3,3n $WORK_DIR/SVC_debug/allDiscordants.us.txt \
    >> $WORK_DIR/SVC_debug/allDiscordants.txt
echo DoneSorting
#Line below is for de novo insertions only
#sort -k 2,2 -k 3,3n $WORK_DIR/allDiscordants.up.us.txt >> $WORK_DIR/allDiscordants.up.txt

time ($SRC/formPEClusters.py -m $7 -p $8 -d ${15} -v ${18} -s ${19} $WORK_DIR/SVC_debug/ \
    $WORK_DIR/SVC_debug/bamStats.txt $WORK_DIR/SVC_debug/binDist.txt)

