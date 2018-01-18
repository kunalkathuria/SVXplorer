#!/bin/bash

FILE_PATH=~/sims/NA12878/BICseq
num=102

for i in $FILE_PATH/*.bin
do
	echo $i
	FILE_TMP="$i.rdtmp"
	rm $FILE_TMP
	touch $FILE_TMP
	while read line
        do
		tmp="$(echo $i | sed -e 's/\/.*\///g')"
		startp="$(echo $line | cut -d " " -f1)"
		endp="$(echo $line | cut -d " " -f2)"
		#echo $startp $endp
		let "eval=$endp-$startp"
		#echo $eval
		if [ "$eval" -lt "$num" ]
			then
				echo $tmp $line >> $FILE_TMP
		fi
        done < $i
done

cat $FILE_PATH/*.rdtmp > ../results/text/RDBin.txt 

