#!/bin/bash
DIR=$2
ST=$4
#$7==pass necessary but didn't work for lumpy-- put as makefile

cat $1 | awk '$7== "PASS"' | grep -v "0\/0" | grep -v "\.\/\." | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DEL | grep -v Cut | awk '{res=$8;split(res,resArray,";"); split(resArray[6],resT,"="); split(resArray[7],resPE, "="); split(resArray[8],resSR,"="); split(resArray[3],len,"="); if (resT[2] >= 6 && (resPE[2] >=4 || resSR[2] >= 4) && (len[2] < -99 || resSR[2] > 0)) print}' > $DIR/DEL.$3.vcf
cat $1 | awk '$7== "PASS"' | grep -v "0\/0" | grep -v "\.\/\." | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DEL | grep Cut | awk '{res=$8;split(res,resArray,";"); split(resArray[9],resT,"="); split(resArray[10],resPE,"="); split(resArray[11],resSR, "="); if (resT[2] >= 6 && (resPE[2] >=4 || resSR[2] >=4)) print}' >> $DIR/DEL.$3.vcf
python ../vcftoBedpe.py $DIR/DEL.$3.vcf 2

cat $1 | awk '$7== "PASS"' | grep -v "0\/0" | grep -v "\.\/\." | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep DUP | grep -v ISINV > $DIR/DUP.$3.vcf
python ../vcftoBedpe.py $DIR/DUP.$3.vcf 2

cat $1 | awk '$7== "PASS"' | grep -v "0\/0" | grep -v "\.\/\." | grep -v GL00 | grep -v NC_ | grep -v hs3 | grep -v -w MT | grep -v BND | grep INV | grep -v ISINV > $DIR/INV.$3.vcf
python ../vcftoBedpe.py $DIR/INV.$3.vcf 2

#cat $1 | awk '$7 == "PASS"' | grep -v "0\/0" | grep -v "\.\/\." | grep BND > $DIR/BND.$3.vcf
#python ../vcftoBedpe.py $DIR/BND.$3.vcf 2
