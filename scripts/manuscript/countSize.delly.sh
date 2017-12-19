#!/bin/bash
rm delcounts.delly.bed
rm tdcounts.delly.bed
rm invcounts.delly.bed

for i in `seq 3 52`;
do
	wc -l $DELLY/src/sim9.50x.1214/delly.sim9.50x.filt.dels.ms.$i.vcf.bed | cut -d " " -f1 >> delcounts.delly.bed
	wc -l $DELLY/src/sim9.50x.1214/delly.sim9.50x.filt.dups.ms.$i.vcf.bed | cut -d " " -f1 >> tdcounts.delly.bed
	wc -l $DELLY/src/sim9.50x.1214/delly.sim9.50x.filt.invs.ms.$i.vcf.bed | cut -d " " -f1 >> invcounts.delly.bed
done
