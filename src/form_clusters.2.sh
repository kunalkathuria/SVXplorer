set -eux

#time ( python ReadDiscordants.py $1 $2 $3 $4 $5 $6 $7 ${13} ${14} ${15} ${16} ${17} ${18} ${19}) # 9 minutes on 50 x homozygous; last fields are map_thresh,AS_THRESH temporary

#echo Sorting
#sort -T ../results/text/ -k 2,2 -k 3,3n ../results/text/All_Discords_P.txt > ../results/text/All_Discords_P_S.txt
#echo DoneSorting

#Line below is for de novo insertions only
#sort -n -k 1,1 ../results/text/All_Discords_I.txt > ../results/text/All_Discords_I_S.txt

time (python FormClusters.py ../results/text/bam_stats.txt $8 $9 ${10} ${11} ${20} ${25} ${24} ${23} ${21} ${22}) #> ../results/text/fc_time.txt # 50m

#time python WriteClusterMap.py $8 # 1m
cp ../results/text/VariantMap.txt ../results/text/VariantMap_O.txt
