MINCS=$1

echo "Classification Start"
cp ../results/text/VariantMap_O.txt ../results/text/VariantMap.txt # uncomment this and following 3 lines when test done
cat ../results/text/All_Clusters.txt | awk '$2 > '$MINCS'' > ../results/text/All_Clusters_minT.txt
sort -k4,4 -k5,5n ../results/text/All_Clusters_minT.txt > ../results/text/All_Clusters_LS.txt
sort -k8,8 -k9,9n ../results/text/All_Clusters_minT.txt > ../results/text/All_Clusters_RS.txt
time (python classifySVs_v28.py $2 $3 $4)

#time python qualifyDeletions.py $7
#mv ../results/text/All_Variants_DC.txt ../results/text/All_Variants.txt

cp ../results/text/ClassifiedVariantMap.txt ../results/text/VariantMap.txt

if [ $6 -eq 1 ]
then
	python SetCover_mq.py $5 $8 ${10} ${11}
else
	python DisjointSetCover.py
fi

python WriteBed.py ../results/text/DisjSetCover_S.txt $7 ${12} ${13}
cp ../results/text/VariantMap_O.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants.txt ../results/text/All_Variants_O.txt

if [ -d ../results/text/pe_results ]
then
	rm -r ../results/text/pe_results
fi
mkdir ../results/text/pe_results
mv ../results/text/*.bedpe ../results/text/pe_results
