#!/bin/bash
DIR=$2
TAG=$3

cat $1 | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DEL > $DIR/DEL.vcf.$TAG
python ../vcftoBedpe.py $DIR/DEL.vcf.$TAG 2

cat $1 | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DUP:TANDEM | > $DIR/DUP.vcf.$TAG
python ../vcftoBedpe.py $DIR/DUP.vcf.$TAG 2

cat $1 | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep INV | grep -v ISINV > $DIR/INV.vcf.$TAG
python ../vcftoBedpe.py $DIR/INV.vcf.$TAG 2

#cat $1 | grep BND > $DIR/BND.vcf.$TAG
#python ../vcftoBedpe.py $DIR/BND.vcf.$TAG 2
