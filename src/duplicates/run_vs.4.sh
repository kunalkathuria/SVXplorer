#cp ../results/text/All_Discords_P_S.sim9.30X ../results/text/All_Discords_P_S.txt
#cp ../results/text/bam_stats.sim9.30X ../results/text/bam_stats.txt 
#cp  ../results/text/bindist.sim9.30X ../results/text/bindist.txt
#cp  ../results/text/VariantMap_SR.sim9.30X ../results/text/VariantMap_SR.txt 
#cp ../results/text/All_Variants_SR.sim9.30X ../results/text/All_Variants_SR.txt
#./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.30X.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.30X.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -C ../data/bed/ceph.hg19.bed -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 30X
#mkdir ../results/text/sim9.30X.1128.vcf
#mv ../results/text/*.bedpe ../results/text/sim9.30X.1128.vcf
#mv ../results/text/All_SVs* ../results/text/sim9.30X.1128.vcf

cp ../results/text/All_Discords_P_S.sim9.5X ../results/text/All_Discords_P_S.txt
cp ../results/text/bam_stats.sim9.5X ../results/text/bam_stats.txt
cp  ../results/text/bindist.sim9.5X ../results/text/bindist.txt
cp  ../results/text/VariantMap_SR.sim9.5X ../results/text/VariantMap_SR.txt
cp ../results/text/All_Variants_SR.sim9.5X ../results/text/All_Variants_SR.txt
./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.5X.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.5X.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -C ../data/bed/ceph.hg19.bed -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 5X
mkdir ../results/text/sim9.5X.1128.vcf
mv ../results/text/*.bedpe ../results/text/sim9.5X.1128.vcf
mv ../results/text/All_SVs* ../results/text/sim9.5X.1128.vcf
cp ../results/text/sim9.5X.1128.vcf/All_SVs.vcf ~/store/results/svc/sim/All_SVs.sim.5x.vcf

cp ../results/text/All_Discords_P_S.sim9.50X.std35.ilb1 ../results/text/All_Discords_P_S.txt
cp ../results/text/bam_stats.sim9.50X.std35.ilb1 ../results/text/bam_stats.txt
cp  ../results/text/bindist.sim9.50X.std35.ilb1 ../results/text/bindist.txt
cp  ../results/text/VariantMap_SR.sim9.50X.std35.ilb1 ../results/text/VariantMap_SR.txt
cp ../results/text/All_Variants_SR.sim9.50X.std35.ilb1 ../results/text/All_Variants_SR.txt
./varsecer -A 0 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/disc.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.50X.bam  -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/sample.50X.bam -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim9/splitters.bam -C ../data/bed/ceph.hg19.bed -D ../results/text/ignoreCHR.txt -F 0 -m 0 -I /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/reference.map.bed 50X 
mkdir ../results/text/sim9.50X.1128.vcf
mv ../results/text/*.bedpe ../results/text/sim9.50X.1128.vcf
mv ../results/text/All_SVs* ../results/text/sim9.50X.1128.vcf
cp ../results/text/sim9.50X.1128.vcf/All_SVs.vcf ~/store/results/svc/sim/All_SVs.sim.50x.vcf
