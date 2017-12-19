#!/bin/bash
for i in `seq 3 52`;
do

	cat lumpy.sim9.50x.$i.gaps.vcf | grep DEL > lumpy.sim9.50x.$i.dels.gaps.vcf
	python ~/vsec/scripts/vcftoBedpe.py lumpy.sim9.50x.$i.dels.gaps.vcf 0
	cat lumpy.sim9.50x.$i.gaps.vcf | grep DUP > lumpy.sim9.50x.$i.dups.gaps.vcf
        python ~/vsec/scripts/vcftoBedpe.py lumpy.sim9.50x.$i.dups.gaps.vcf 0
	cat lumpy.sim9.50x.$i.gaps.vcf | grep INV > lumpy.sim9.50x.$i.invs.gaps.vcf
        python ~/vsec/scripts/vcftoBedpe.py lumpy.sim9.50x.$i.invs.gaps.vcf 0
done
