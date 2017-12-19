#!/bin/bash
MINCS=$1
WORK_DIR=${15}

echo "Classification Start"
cp $WORK_DIR/VariantMap_O.txt $WORK_DIR/VariantMap.txt # uncomment this and following 3 lines when test done
cat $WORK_DIR/All_Clusters.txt | awk '$2 >= '$MINCS'' > $WORK_DIR/All_Clusters_minT.txt
sort -k4,4 -k5,5n $WORK_DIR/All_Clusters_minT.txt > $WORK_DIR/All_Clusters_LS.txt
sort -k7,7 -k8,8n $WORK_DIR/All_Clusters_minT.txt > $WORK_DIR/All_Clusters_RS.txt
time (python classifySVs_v28.py $2 $3 $4 $WORK_DIR)

#time python qualifyDeletions.py $7
#mv $WORK_DIR/All_Variants_DC.txt $WORK_DIR/All_Variants.txt

cp $WORK_DIR/ClassifiedVariantMap.txt $WORK_DIR/VariantMap.txt

if [ $6 -eq 1 ]
then
	python SetCover_mq.py $5 $8 ${10} ${11} ${14} $WORK_DIR
else
	python DisjointSetCover.py $WORK_DIR
fi

python WriteBed.o2.py $WORK_DIR/DisjSetCover_S.txt $WORK_DIR
cp $WORK_DIR/VariantMap_O.txt $WORK_DIR/VariantMap.txt
cp $WORK_DIR/All_Variants.txt $WORK_DIR/All_Variants_O.txt

if [ -d $WORK_DIR/pe_results ]
then
	rm -r $WORK_DIR/pe_results
fi
mkdir $WORK_DIR/pe_results
mv $WORK_DIR/*.bedpe $WORK_DIR/pe_results
