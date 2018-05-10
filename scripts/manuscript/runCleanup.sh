#!/usr/bin/env bash
wdir=$3
clPath=$2
ovBuffer=200
igBuffer=100
mergeBuffer=100
sampleBAM=$1
suppPerc=.5
suppMin=5

sort -k4,4 -k5,5n $clPath > ${clPath}.ls
clPathN=${clPath}.ls
python scripts/manuscript/pickCleanClusters.str.sprop.py $clPathN $wdir $ovBuffer
sort -k4,4 -k5,5n $wdir/allClusters.clean.str.txt > $wdir/allClusters.clean.str.ls2.txt
cp  $wdir/badRegions.bed  $wdir/br1.bed
python scripts/manuscript/pickCleanClusters.allO.sprop.py $wdir/allClusters.clean.str.ls2.txt $wdir
mv  $wdir/badRegions.bed  $wdir/br2.bed
cat  $wdir/br1.bed  $wdir/br2.bed >  $wdir/badRegions.bed
sort -k1,1 -k2,2n $wdir/badRegions.bed > $wdir/badRegions.sorted.bed
bedtools merge -d $mergeBuffer -i $wdir/badRegions.sorted.bed > $wdir/badRegions.m.bed
sort -k4,4 -k5,5n $wdir/allClusters.clean.allO.txt > $wdir/allClusters.clean.allO.ls.txt
python scripts/manuscript/postCleanup.sprop.py $wdir/allClusters.clean.allO.ls.txt $wdir/badRegions.m.bed \
    $wdir $sampleBAM $igBuffer $suppPerc $suppMin
echo Cleaned up
#rm $clPathN $wdir/br1.bed $wdir/br2.bed $wdir/badRegions.sorted.bed \
#    $wdir/allClusters*clean*txt
