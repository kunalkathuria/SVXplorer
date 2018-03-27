#!/bin/bash
DIR=$2
TAG=$3
cat $1 | grep DEL > $DIR/DEL.vcf.$TAG
python ../vcftoBedpe.py $DIR/DEL.vcf.$TAG 0

cat $1 | grep DUP > $DIR/DUP.vcf.$TAG
python ../vcftoBedpe.py $DIR/DUP.vcf.$TAG 0

cat $1 | grep INV > $DIR/INV.vcf.$TAG
python ../vcftoBedpe.py $DIR/INV.vcf.$TAG 0

cat $1 | grep BND > $DIR/BND.vcf.$TAG
python ../vcftoBedpe.py $DIR/BND.vcf.$TAG 0
