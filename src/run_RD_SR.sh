#samtools view -h ${14} \
#    | ~/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
#    | samtools view -Sb - \
#    > /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/splitters.us.bam

#samtools sort -n -@32 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/splitters.bam
#rm /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/splitters.us.bam

#cp ../results/text/ClassifiedVariantMap.txt ../results/text/VariantMap.txt
#cp ../results/text/All_Variants_O.txt ../results/text/All_Variants.txt
#python add_SR_hash.2.py $1 $2 $3 $4 ${13} ${12} ${10} ${18} ${21}
#cp ../results/text/VariantMap_SR.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_SR.txt ../results/text/All_Variants.txt
#python SetCover_mq.py $6 ${11} ${15} ${16}
#python WriteBed.o2.py ../results/text/DisjSetCover_S.txt
#if [ -d ../results/text/sr_results ]
#then
#        rm -r ../results/text/sr_results
#fi
#mkdir ../results/text/sr_results
#mv ../results/text/*.bedpe ../results/text/sr_results

if [ $7 -eq 1 ]
then
	counter=0
	while read -r line;
	do
        	counter=$((counter+1))
        	if [[ $counter -eq 2 ]]; then
                	echo $line
                	meanIL=${line%.*}
        	elif [[ $counter -eq 6 ]]; then
                	echo ${line%.*}
                	DISC_D=${line%.*}
        	fi
	done < <(tr -d '\r' < ../results/text/bam_stats.txt)
	RD_SLOP=$((meanIL+DISC_D ))
	echo $RD_SLOP
	cat ../results/text/All_Variants.txt | grep -v INV | grep -v known > ../results/text/inDels.txt # all INDELS
	sort -k6,6 ../results/text/inDels.txt > ../results/text/inDels_S.txt
	#./formatRD.sh #make first 2 vars in here user-input
	python add_RD.o2.py $5 $8 ${15} ${16} ${17} ${19} $RD_SLOP ${24} ${23}
	cp ../results/text/VariantMap_RD.txt ../results/text/VariantMap.txt
	cp ../results/text/All_Variants_RD.txt ../results/text/All_Variants.txt
fi

#if [ $9 -eq 1 ]
#then
#	python SetCover_mq.py $6 ${11} ${15} ${16}
#else
#	python DisjointSetCover.py
#fi

python WriteBed.o5.py ../results/text/DisjSetCover_S.txt ${14} ${22} ${23} ${25} ${26} ${27} ${28} ${29}
cat ../results/text/deletions.bedpe ../results/text/tandemDuplications.bedpe ../results/text/inversions.bedpe ../results/text/insertions.bedpe ../results/text/unknowns.bedpe > ../results/text/All_SVs
python ../scripts/bed2vcf.py ../results/text/All_SVs
cp ../results/text/VariantMap_O.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_O.txt ../results/text/All_Variants.txt

