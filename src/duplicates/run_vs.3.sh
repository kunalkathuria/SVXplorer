./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.10X.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.10X.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 10X
mkdir ../results/text/sim9.10X
mv ../results/text/*.bedpe ../results/text/sim9.10X
mv ../results/text/All_SVs* ../results/text/sim9.10X
mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim9.10X
mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim9.10X
mv ../results/text/bindist.txt ../results/text/bindist.sim9.10X
mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim9.10X
mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim9.10X

