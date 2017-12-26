#!/bin/bash
set -eux
WORK_DIR=${26}
SRC=${27}

time ( python $SRC/ReadDiscordants.graph.py $1 $2 $3 $4 $5 $6 $7 ${13} ${14} ${15} ${16} ${17} ${18} ${19} $WORK_DIR/SVC_debug) # 9 minutes on 50 x homozygous; last fields are map_thresh,AS_THRESH temporary
echo Sorting
sort -T $WORK_DIR/SVC_debug -k 2,2 -k 4,4 -k 3,3n $WORK_DIR//SVC_debug/allDiscordants.us.txt > $WORK_DIR/SVC_debug/allDiscordants.txt
echo DoneSorting

#Line below is for de novo insertions only
#sort -n -k 1,1 $WORK_DIR/allDiscordants.up.us.txt > $WORK_DIR/allDiscordants.up.txt

time (python $SRC/FormClusters.graph.sec.py $WORK_DIR/SVC_debug/bamStats.txt $8 $9 ${10} ${11} ${20} ${25} ${24} ${23} ${21} ${22} $WORK_DIR/SVC_debug) #> $WORK_DIR/fc_time.txt # 50m

