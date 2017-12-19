#!/bin/bash
rm delcounts.manta.bed
rm tdcounts.manta.bed
rm invcounts.manta.bed

for i in `seq 3 52`;
do
	wc -l ~/store/results/manta/sim9.50x/results/variants/DEL.$i.vcf.bed | cut -d " " -f1 >> delcounts.manta.bed
	wc -l ~/store/results/manta/sim9.50x/results/variants/DUP.$i.vcf.bed| cut -d " " -f1 >> tdcounts.manta.bed
	wc -l ~/store/results/manta/sim9.50x/results/variants/INV.$i.vcf.bed | cut -d " " -f1 >> invcounts.manta.bed
done
