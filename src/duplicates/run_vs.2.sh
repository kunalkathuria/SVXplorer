#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.1X.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.1X.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -C ../data/bed/ceph.hg19.bed -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 1X
#mkdir ../results/text/sim9.1X
#mv ../results/text/*.bedpe ../results/text/sim9.1X
#mv ../results/text/All_SVs* ../results/text/sim9.1X
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim9.1X
#mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim9.1X
#mv ../results/text/bindist.txt ../results/text/bindist.sim9.1X
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim9.1X
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim9.1X

./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -C ../data/bed/ceph.hg19.bed -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 50X
mkdir ../results/text/sim9.50X.std35.ilb1
mv ../results/text/*.bedpe ../results/text/sim9.50X.std35
mv ../results/text/All_SVs* ../results/text/sim9.50X.std35
mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.sim9.50X.std35.ilb1
mv ../results/text/bam_stats.txt ../results/text/bam_stats.sim9.50X.std35.ilb1
mv ../results/text/bindist.txt ../results/text/bindist.sim9.50X.std35.ilb1
mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.sim9.50X.std35.ilb1
mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.sim9.50X.std35.ilb1
