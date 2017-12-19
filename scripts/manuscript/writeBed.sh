#!/bin/bash
for i in `seq 3 52`;
do

	cat ./results/variants/diploidSV.vcf.minsupp.$i | grep DEL > ./results/variants/DEL.$i.vcf
	python ~/vsec/scripts/vcftoBedpe.py ./results/variants/DEL.$i.vcf 0
	cat ./results/variants/diploidSV.vcf.minsupp.$i | grep DUP > ./results/variants/DUP.$i.vcf
	python ~/vsec/scripts/vcftoBedpe.py ./results/variants/DUP.$i.vcf 0
	cat ./results/variants/diploidSV.vcf.minsupp.$i | grep INV > ./results/variants/INV.$i.vcf
	python ~/vsec/scripts/vcftoBedpe.py ./results/variants/INV.$i.vcf 0

done
