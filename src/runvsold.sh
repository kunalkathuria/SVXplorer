#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.25.cov.15.bam  -b ~/vsec/data/bams/sim2.stdev.25.cov.15.bam -r splitters.25.15.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.R1
#mv ../results/text/sr_results ../results/text/sr_results.R1
#mkdir ../results/text/R1
#mv ../results/text/*.bedpe ../results/text/R1
#
#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.25.cov.50.sorted.bam  -b ~/vsec/data/bams/sim2.stdev.25.cov.50.bam -r splitters.25.50.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.R2
#mv ../results/text/sr_results ../results/text/sr_results.R2
#mkdir ../results/text/R2
#mv ../results/text/*.bedpe ../results/text/R2
#
#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.25.cov.70.sorted.bam  -b ~/vsec/data/bams/sim2.stdev.25.cov.70.sorted.bam -r splitters.25.70.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.R3
#mv ../results/text/sr_results ../results/text/sr_results.R3
#mkdir ../results/text/R3
#mv ../results/text/*.bedpe ../results/text/R3
#
#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.50.cov.50.sorted.bam  -b ~/vsec/data/bams/sim2.stdev.50.cov.50.sorted.bam -r splitters.50.50.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.R4
#mv ../results/text/sr_results ../results/text/sr_results.R4
#mkdir ../results/text/R4
#mv ../results/text/*.bedpe ../results/text/R4
#
#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.5.cov.50.sorted.bam  -b ~/vsec/data/bams/sim2.stdev.5.cov.50.sorted.bam -r splitters.5.50.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.R5
#mv ../results/text/sr_results ../results/text/sr_results.R5
#mkdir ../results/text/R5
#mv ../results/text/*.bedpe ../results/text/R5

#./varsecer -A 0 -B 32 -z  /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/AJtrio/HG002/alignments/discordant.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/AJtrio/HG002/alignments/sample.bam -b /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/AJtrio/HG002/alignments/sample.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/AJ/HG002/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed

#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/srr505.final.cnv -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.srr505.prim.wt.t017
#mv ../results/text/sr_results ../results/text/sr_results.srr505.prim.wt.t017
#mkdir ../results/text/srr505.prim.wt.t017
#mv ../results/text/*bedpe* ../results/text/srr505.prim.wt.t017

./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sec.err
#mv ../results/text/sr_results ../results/text/sr_results.s9.sec.err
#mkdir ../results/text/s9.sec.err
#mv ../results/text/*bedpe* ../results/text/s9.sec.err
#mv ../results/text/All_SVs* ../results/text/s9.sec.err
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.err.sec.txt
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.err.txt

#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/disc.sim100.all.bam -a ~/vsec/data/bams/sim2.stdev.25.cov.70.sorted.bam  -b ~/vsec/data/bams/sim2.stdev.25.cov.70.sorted.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim2/sim100.all.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim2.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sim2.prim
#mv ../results/text/sr_results ../results/text/sr_results.s9.sim2.prim
#mkdir ../results/text/s9.sim2.prim
#mv ../results/text/*bedpe* ../results/text/s9.sim2.prim
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim2.prim.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim2.s9.prim.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim2.s9.prim.txt
##mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim2.txt
#
#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/disc.sim100.all.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim4.final.cnv -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sim4.prim
#mv ../results/text/sr_results ../results/text/sr_results.s9.sim4.prim
#mkdir ../results/text/s9.sim4.prim
#mv ../results/text/*bedpe* ../results/text/s9.sim4.prim
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim4.prim.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim4.s9.prim.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim4.s9.prim.txt
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim4.txt

#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.err.cc.1
#mv ../results/text/sr_results ../results/text/sr_results.err.cc.1
#mkdir ../results/text/err.cc.1
#mv ../results/text/*bedpe* ../results/text/err.cc.1
#mv ../results/text/All_SVs* ../results/text/err.cc.1
#cp ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.err.prim.txt
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.err.txt

#./varsecer2 -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed

#./varsecer -A 0 -B 32 -a ~/vsec/data/bams/del.800bp.test.sorted.bam -b ~/vsec/data/bams/del.800bp.test.sorted.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -F 0 -m 0

#./varsecer -A 0 -B 32 -z  /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/discordant.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed0
