#!/bin/bash
DIR=$2

cat $1 | awk '$7 == "PASS"' | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DEL > $DIR/DEL.vcf
python ../vcftoBedpe.py $DIR/DEL.vcf 2

cat $1 | awk '$7 == "PASS"' | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DUP:TANDEM > $DIR/DUP.vcf
python ../vcftoBedpe.py $DIR/DUP.vcf 2

cat $1 | awk '$7 == "PASS"' | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep INV | grep -v ISINV > $DIR/INV.vcf
python ../vcftoBedpe.py $DIR/INV.vcf 2

#cat $1 | awk '$7 == "PASS"' | grep BND > $DIR/BND.vcf
#python ../vcftoBedpe.py $DIR/BND.vcf 2
