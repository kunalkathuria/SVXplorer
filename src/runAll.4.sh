#cp ../results/text/All_Discords_P_S.s2.ERR194147.discd95.txt ../results/text/All_Discords_P_S.txt
#./run_sv_caller.14.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt
#mkdir ../results/text/s2.srr505.19.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.19.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.19
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.19.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.19.txt

#cp ../results/text/All_Discords_P_S.s2.sim4.discd9985.txt ../results/text/All_Discords_P_S.txt
#cp ../results/text/All_Variants_SR.sim4_s2.txt ../results/text/All_Variants_SR.txt
#cp ../results/text/VariantMap_SR.sim4_s2.txt ../results/text/VariantMap_SR.txt
#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim4_3.nogaps.cnv
#mkdir ../results/text/s2.srr505.20.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.20.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.20
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.20.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.20.txt

#./run_cnv.filters.sh sim4
#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/disc.sim100.all.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim4.nogaps.cnv
#mkdir ../results/text/s2.sim4
#mv ../results/text/*.bedpe ../results/text/s2.sim4
#mv ../results/text/pe_results ../results/text/pe_results.s2.sim4
#cp ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.sim4.txt
#cp ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.sim4.txt
#cp ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.s2.sim4.discd9985.txt

./run_cnv.filters.sh sim4 /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam hg19 1

./run_cnv.filters.sh srr505 /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam GRCh37 2
cp ../results/text/All_Discords_P_S.s2.srr505.discd9985.txt ../results/text/All_Discords_P_S.txt
./run_sv_caller.2.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/srr505.final.cnv

#mkdir ../results/text/s2.srr505.15.rd
mv ../results/text/*.bedpe ../results/text/s2.srr505.15.rd
mv ../results/text/pe_results/* ../results/text/pe_results.s2.srr505.15
mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.15.txt
mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.15.txt

./run_cnv.filters.sh srr339 /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR3397076/alignments/sample.bam GRCh37 2
./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR3397076/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/srr339.final.cnv

mkdir ../results/text/s2.srr339.18.rd #18 means .60 dd here
mv ../results/text/*.bedpe ../results/text/s2.srr339.18.rd
mv ../results/text/pe_results ../results/text/pe_results.s2.srr339.18

./run_cnv.filters.sh err /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam GRCh37 2
#cp ../results/text/All_Discords_P_S.s2.err.discd95.txt ../results/text/All_Discords_P_S.txt
#./run_sv_caller.2.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv

#mkdir ../results/text/s2.err.18.rd
#mv ../results/text/*.bedpe ../results/text/s2.err.18.rd
#mv ../results/text/pe_results/* ../results/text/pe_results.s2.err.18

#./run_sv_caller.14.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt

#mkdir ../results/text/s2.srr505.4.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.4.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.4
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.4.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.4.txt

#./run_sv_caller.2.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt

#mkdir ../results/text/s2.srr505.18.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.18.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.18
#cp ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.s2.ERR194147.discd95.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.18.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.18.txt

#./run_sv_caller.3.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt

#mkdir ../results/text/s2.srr505.8.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.8.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.8
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.8.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.8.txt

#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt

#mkdir ../results/text/s2.srr505.15.nord
#mv ../results/text/*.bedpe ../results/text/s2.srr505.15.nord
#mv ../results/text/pe_results ../results/text/pe_results.s2.srr505.15
#cp ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.s2.ERR194147.discd9985.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.s2.srr505.15.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.s2.srr505.15.txt
