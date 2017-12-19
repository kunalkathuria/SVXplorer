#!/bin/bash
for i in `seq 3 52`;
do

	cat sim9.50x.ms.$i/deletions.bedpe | cut -f1,2,6 > sim9.50x.ms.$i/deletions.bed
	cat sim9.50x.ms.$i/tandemDuplications.bedpe | cut -f1,2,6 > sim9.50x.ms.$i/tandemDuplications.bed
	cat sim9.50x.ms.$i/inversions.bedpe | cut -f1,2,6 > sim9.50x.ms.$i/inversions.bed
done
