REMOTE_FILE1=/m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/homozygous_variants/alignments/id.ns.alignment.bam
REMOTE_FILE2=/m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/homozygous_variants/alignments/id.alignment.bam

#echo $1 $2 $3
time python ReadDiscordants.py $REMOTE_FILE1 $REMOTE_FILE2 20 1 .95 500000 0 # 9 minutes on 50 x homozygous

echo Sorting
sort -T ../results/text/ -k 2,2 -k 3,3n ../results/text/All_Discords_P.txt > ../results/text/All_Discords_P_S.txt
echo DoneSorting

#Line below is for random insertions only
#sort -n -k 1,1 ../results/text/All_Discords_I.txt > ../results/text/All_Discords_I_S.txt

time (python FormClusters.py ../results/text/bam_stats.txt 4 5 20 2) # 50m
time (sort -T ../results/text/ -n -k 2,2 ../results/text/VariantMapInp_P.txt > ../results/text/VariantMapInp.txt) #1m
time python WriteClusterMap.py 4 # 1m
cp ../results/text/VariantMap.txt ../results/text/VariantMap_O.txt
