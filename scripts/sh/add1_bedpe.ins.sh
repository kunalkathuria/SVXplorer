file=$1
cat $file | cut -f1 > c1
cat $file | cut -f5 > c2
cat $file | cut -f6 > c6
cat $file | awk '{print $5+1}' > c3
cat $file | awk '{print $6-1}' > c5
paste c1 c2 c3 c1 c5 c6 > $file.6col.bedpe
rm c1 c2 c3 c5 c6
