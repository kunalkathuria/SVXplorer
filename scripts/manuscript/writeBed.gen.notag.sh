#!/bin/bash
DIR=$2

cat $1 | grep PASS | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DEL > $DIR/DEL.2.vcf
python ../vcftoBedpe.py $DIR/DEL.2.vcf 2

cat $1 | grep PASS | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DUP:TANDEM > $DIR/DUP.2.vcf
python ../vcftoBedpe.py $DIR/DUP.2.vcf 2

cat $1 | grep PASS | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep INV | grep -v ISINV > $DIR/INV.2.vcf
python ../vcftoBedpe.py $DIR/INV.2.vcf 2

#cat $1 | grep PASS | grep BND > $DIR/BND.2.vcf
#python ../vcftoBedpe.py $DIR/BND.2.vcf 2
