#!/bin/bash
rm delcounts.bed
rm tdcounts.bed
rm invcounts.bed

for i in `seq 3 52`;
do
	wc -l ~/vsec/results/text/sim9.50x.ms.$i/deletions.bed | cut -d " " -f1 >> delcounts.bed
	wc -l ~/vsec/results/text/sim9.50x.ms.$i/tandemDuplications.bed | cut -d " " -f1 >> tdcounts.bed
	wc -l ~/vsec/results/text/sim9.50x.ms.$i/inversions.bed | cut -d " " -f1 >> invcounts.bed
done
