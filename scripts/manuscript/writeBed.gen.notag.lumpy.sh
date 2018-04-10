#!/bin/bash
DIR=$2
#$7==pass necessary but didn't work for lumpy-- put as makefile

cat $1 | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep -v "0\/0" | grep -v "\.\/\." | grep DEL > $DIR/DEL.$3.vcf
python ../vcftoBedpe.py $DIR/DEL.$3.vcf 2

cat $1 | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep -v "0\/0" | grep -v "\.\/\." | grep DUP > $DIR/DUP.$3.vcf
python ../vcftoBedpe.py $DIR/DUP.$3.vcf 2

cat $1 | grep -v "\.\/\." | grep -v "0\/0" | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep INV | grep -v ISINV > $DIR/INV.$3.vcf
python ../vcftoBedpe.py $DIR/INV.$3.vcf 2

#cat $1 | awk '$7 == "PASS"' | grep BND > $DIR/BND.$3.vcf
#python ../vcftoBedpe.py $DIR/BND.$3.vcf 2
