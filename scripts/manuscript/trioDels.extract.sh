cat deletions.bedpe | grep "1/1" > HZ.dels.bedpe
cat ~/vscopy/results/text/$1/deletions.bedpe ~/vscopy/results/text/$1/deletions.bedpe > AJ.bothparents.dels.bedpe
cat HZ.dels.bedpe | cut -f1-6 > HZ.dels.2.bedpe
cat AJ.bothparents.dels.bedpe | cut -f1-6 > AJ.bothparents.dels.2.bedpe
bedtools pairtopair -type both -is -slop 0 -a HZ.dels.2.bedpe -b AJ.bothparents.dels.2.bedpe | cut -f1-6 | sort | uniq |wc -l
wc -l HZ.dels.2.bedpe
