#!/bin/bash
cp ../results/text/VariantMap_SR.txt ../results/text/VariantMap.txt
cp ../results/text/All_Variants_SR.txt ../results/text/All_Variants.txt

for i in `seq 3 53`;
do

	if [ $9 -eq 1 ]
	then
		python SetCover_mq.loop.py $6 ${11} ${15} ${16} ${31} $i $i $i
	else
		python DisjointSetCover.py
	fi

	python WriteBed.o8.o.py ../results/text/DisjSetCover_S.txt ${14} ${22} ${23} ${25} ${26} ${27} ${28} ${29} ${30}
	echo "WB Done"
	cat ../results/text/deletions.bedpe ../results/text/tandemDuplications.bedpe ../results/text/inversions.bedpe ../results/text/insertions.bedpe ../results/text/unknowns.bedpe > ../results/text/All_SVs
	python ../scripts/bed2vcf.py ../results/text/All_SVs
	mkdir ../results/text/sim8.50x.ms.$i
	mv ../results/text/*.bedpe ../results/text/sim8.50x.ms.$i
	mv ../results/text/All_SVs* ../results/text/sim8.50x.ms.$i
done
