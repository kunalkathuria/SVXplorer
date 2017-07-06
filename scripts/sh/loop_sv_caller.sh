set -eux

./run_sv_caller.2.sh -A 1 -B 32 -z /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.discordant.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.alignment.bam -b /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.alignment.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt
mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.mq1ras95.txt
mv ../results/text/pe_results ../results/text/pe_results_mq1
mkdir ../results/text/res_mq1
mv ../results/text/*.bedpe ../results/text/res_mq1
./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.discordant.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.alignment.bam -b /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.alignment.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt
mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.mq0AS0.txt
mv ../results/text/pe_results ../results/text/pe_results_mq0_aas0
mkdir ../results/text/res_mq0_aas0
mv ../results/text/*.bedpe ../results/text/res_mq0_aas0
#./run_sv_caller.2.sh -A 1 -B 32 -z /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.discordant.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.alignment.bam -b /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.alignment.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt
#mv ../results/text/pe_results ../results/text/pe_results_f
#mkdir ../results/text/res_sec_mq0_95
#mv ../results/text/*.bedpe ../results/text/res_sec_mq0_95

#res1 20 .85 20 (AS_REL = .95)
#res2 20 .95 20 (AS_REL = .95)
#res3 20 0 20 (AS_REL = 0)
#res4 mq 20 .85 abs ms .95 rel ms
#res5 mq 0 .85 abs as .95 rel as
#res6 mq 0 .85 abs ms .95 rel ms
#res7 mq 0 .95 as .95 rel as
#res8 mq 0 .95 ms .95relms
