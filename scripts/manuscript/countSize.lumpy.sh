#!/bin/bash
rm delcounts.lumpy.bed
rm tdcounts.lumpy.bed
rm invcounts.lumpy.bed

for i in `seq 3 52`;
do
	wc -l /home/kk7t/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/lumpy.sim9.50x.$i.dels.gaps.vcf.bed | cut -d " " -f1 >> delcounts.lumpy.bed
	wc -l /home/kk7t/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/lumpy.sim9.50x.$i.dups.gaps.vcf.bed | cut -d " " -f1 >> tdcounts.lumpy.bed
	wc -l /home/kk7t/store/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/lumpy.sim9.50x.$i.invs.gaps.vcf.bed| cut -d " " -f1 >> invcounts.lumpy.bed
done
