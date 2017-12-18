#!/bin/bash
set -eux
WORK_DIR=${26}

time ( python ReadDiscordants.graph.py $1 $2 $3 $4 $5 $6 $7 ${13} ${14} ${15} ${16} ${17} ${18} ${19} $WORK_DIR) # 9 minutes on 50 x homozygous; last fields are map_thresh,AS_THRESH temporary
echo Sorting
sort -T $WORK_DIR -k 2,2 -k 4,4 -k 3,3n $WORK_DIR/All_Discords_P.txt > $WORK_DIR/All_Discords_P_S.txt
echo DoneSorting

#Line below is for de novo insertions only
#sort -n -k 1,1 $WORK_DIR/All_Discords_I.txt > $WORK_DIR/All_Discords_I_S.txt

time (python FormClusters.graph.sec.py $WORK_DIR/bam_stats.txt $8 $9 ${10} ${11} ${20} ${25} ${24} ${23} ${21} ${22} $WORK_DIR) #> $WORK_DIR/fc_time.txt # 50m

#time python WriteClusterMap.py $8 # 1m
cp $WORK_DIR/VariantMap.txt $WORK_DIR/VariantMap_O.txt
