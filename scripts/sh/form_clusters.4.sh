set -eux

#python ~/store/kk7t/sv_caller/gitSVCaller/code/ReadDiscordants.py $1 $2 20 1 .95 500000 0
#time ( python ReadDiscordants.py $1 $2 $3 $4 $5 $6 $7 ${13} ${14} ${15} ${16} ${17} ${18} .999) # 9 minutes on 50 x homozygous; last fields are map_thresh,AS_THRESH temporary

#echo Sorting
#sort -T ../results/text/ -k 2,2 -k 3,3n ../results/text/All_Discords_P.txt > ../results/text/All_Discords_P_S.txt
#echo DoneSorting

#Line below is for random insertions only
#sort -n -k 1,1 ../results/text/All_Discords_I.txt > ../results/text/All_Discords_I_S.txt

time (python FormClusters.py ../results/text/bam_stats.txt $8 $9 ${10} ${11} .9999 1.0 1 1.67 20 2.0) #> ../results/text/fc_time.txt # 50m
#time (sort -T ../results/text/ -n -k 2,2 ../results/text/VariantMapInp_P.txt > ../results/text/VariantMapInp.txt) #1m

#time python WriteClusterMap.py $8 # 1m
cp ../results/text/VariantMap.txt ../results/text/VariantMap_O.txt
#python SetCover_mq.py $8 ${17} 3000000 100000000
