#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/srr505.final.cnv -F 0 
#mv ../results/text/pe_results ../results/text/pe_results.s4.srr505.15
#mv ../results/text/sr_results ../results/text/sr_results.s4.srr505.15 
#mkdir ../results/text/s6.srr505.15
#mv ../results/text/*bedpe* ../results/text/s6.srr505.15

#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sec.err
#mv ../results/text/sr_results ../results/text/sr_results.s9.sec.err
#mkdir ../results/text/s9.sec.err
#mv ../results/text/*bedpe* ../results/text/s9.sec.err
#mv ../results/text/All_SVs* ../results/text/s9.sec.err
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.err.sec.txt
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.err.txt

#./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim3/disc.sim100.all.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim3/sim100.all.sorted.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim3/sim100.all.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim3/sim100.all.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/sim3.final.cnv  1.0 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sim3.prim
#mv ../results/text/sr_results ../results/text/sr_results.s9.sim3.prim
#mkdir ../results/text/s9.sim3.prim
#mv ../results/text/*bedpe* ../results/text/s9.sim3.prim
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim3.prim.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim3.s9.prim.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim3.s9.prim.txt
##mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim3.txt
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

./varsecer -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR194147/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/err.final.cnv -G .95 -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed
#mv ../results/text/pe_results ../results/text/pe_results.s9.sec.err
#mv ../results/text/sr_results ../results/text/sr_results.s9.sec.err
#mkdir ../results/text/s9.sec.err
#mv ../results/text/*bedpe* ../results/text/s9.sec.err
#mv ../results/text/All_SVs* ../results/text/s9.sec.err
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.err.sec.txt
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.err.txt
