set -eux

BWA=~/scratch/share/bwa-0.7.12/bwa
#SAMTOOLS=~/samtools-1.3/samtools
SAMTOOLS=/m/cphg-RLscratch/cphg-RLscratch/share/samtools-0.1.19/samtools
WGSIM=~/scratch/share/wgsim/wgsim
SVSIM=./sv_sim/SVsim
SVCommands=./sv_sim/SVcommands_test_1SV.sim
REFERENCE=/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37_decoy.fasta
#REFERENCE=/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta
#REFERENCE=../data/refs/ref_1k.fa
TARGET=/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim5/genome_rearranged.fasta #../data/refs/targ_del.fa
TARGETNAME=../data/sims/svsim_out/test_1SV.root
READ1=/home/kk7t/store/../target/NA12878/SRR3397076/1_ed.fq #/m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/data/ERR194147_pass_1.fastq
READ2=/home/kk7t/store/../target/NA12878/SRR3397076/2_ed.fq #/m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/data/ERR194147_pass_2.fastq
 #/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim5/sim100.all.1.fq
DATA=../other_tools/lumpy/data
REMOTE_FILE=~/homozygous/alignments/id.ns.alignment.bam

MEAN=350
STDEV=50
RDL=100
NBP=3100000000
COVERAGE=50
NREADS=$(( (COVERAGE * NBP) / (2 * RDL) ))
NSIMS=1
MARGIN=200

echo $NREADS reads
# prepare simualted target
#python $SVSIM -i $SVCommands -r $REFERENCE -o $TARGETNAME -W -d -n $NSIMS

# prep fake reads from the target
#$WGSIM -N$NREADS -1 $RDL -2 $RDL -d $MEAN -s $STDEV $TARGET $READ1 $READ2

# bwa alignment
#$BWA index $REFERENCE
#python convertFQIllumina.py /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/data/SRR505885_pass_1.fastq.gz /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/data/SRR505885_pass_2.fastq.gz  /home/kk7t/store/../target/NA12878/SRR505885/
#$BWA mem -R '@RG\tID:foo\tSM:bar' -a -Y -t 32 $REFERENCE $READ1 $READ2 \
#| $SAMTOOLS view -S -b - \
#> /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/nsall.bam
#samtools sort -@32 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim5/sim100.all.sorted.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim5/sim100.all.bam
#samtools index /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim5/sim100.all.sorted.bam

#mv ../results/text/VariantMap_SR.srr339.9.txt ../results/text/VariantMap_SR.txt
#mv ../results/text/All_Variants_SR.srr339.9.txt ../results/text/All_Variants_SR.txt
#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR3397076/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr339.txt
#mkdir ../results/text/sr.srr339.9.rd
#mv ../results/text/*.bedpe ../results/text/sr.srr339.9.rd
#mv ../results/text/pe_results ../results/text/pe_results.srr339.9.rd
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.srr339.9.rd.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.srr339.9.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.srr339.9.txt

#mv ../results/text/All_Discords_P_S.srr505.9.txt ../results/text/All_Discords_P_S.txt
#mv ../results/text/VariantMap_SR.srr505.9.txt ../results/text/VariantMap_SR.txt
#mv ../results/text/All_Variants_SR.srr505.9.txt ../results/text/All_Variants_SR.txt
./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/SRR505885/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt -x ../results/text/cnvs.srr505.txt
#mkdir ../results/text/sr.srr505.9_sronly
#mv ../results/text/*.bedpe ../results/text/sr.srr505.9_sronly
#mv ../results/text/pe_results ../results/text/pe_results.srr505.9.rd
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.srr505.9.rd.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.srr505.9.rd.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.srr505.9.rd.txt
#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR/disc.all.bam -a /m/cphg-RLscratch/cphg-RLscratch/ar7jq/read_depth/NA12878/ERR194147/alignments/sample.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR/nsall.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/ERR/splitters.ns.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt
#mkdir ../results/text/sr.err.9.prim
#mv ../results/text/*.bedpe ../results/text/sr.err.9.prim
#mv ../results/text/All_Discords_P_S.txt ../results/text/All_Discords_P_S.err.9.prim.txt
#mv ../results/text/VariantMap_SR.txt ../results/text/VariantMap_SR.err.9.prim.txt
#mv ../results/text/All_Variants_SR.txt ../results/text/All_Variants_SR.err.9.prim.txt
#time ($SAMTOOLS view -F 3586 -b -o ../data/bams/disc.bam ../data/bams/aln_test.bam)
#time ($SAMTOOLS view -F 3586 -b -o ../data/bams/disc.bam $REMOTE_FILE) #34m
#time ($SAMTOOLS view -f 64 -b -o ../data/bams/aln1s_100.bam ../data/bams/disc.bam) #4m
#time ($SAMTOOLS view -f 128 -b -o ../data/bams/aln2s_100.bam ../data/bams/disc.bam) # 4m

# run set cover
#./runSetCover7_indiv_temp.sh
#python formatBP_chrNum_v2.py ${TARGETNAME}.bedpe
#./run_classify.sh 25
#./run_PE_RD.sh $7 $8
