set -eux

BWA=~/scratch/share/bwa-0.7.12/bwa
#SAMTOOLS=~/samtools-1.3/samtools
SAMTOOLS=/m/cphg-RLscratch/cphg-RLscratch/share/samtools-0.1.19/samtools
WGSIM=~/scratch/share/wgsim/wgsim
SVSIM=./sv_sim/SVsim
SVCommands=./sv_sim/SVcommands_test_1SV.sim
REFERENCE=/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta
#REFERENCE=../data/refs/ref_1k.fa
TARGET=/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/genome_rearranged.fasta #../data/refs/targ_del.fa
TARGETNAME=../data/sims/svsim_out/test_1SV.root
READ1=/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.1.fq
READ2=/m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.2.fq
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
$BWA mem -R '@RG\tID:foo\tSM:bar' -a -Y -t 19 $REFERENCE $READ1 $READ2 \
| $SAMTOOLS view -S -b - \
> /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.bam
samtools sort -@18 -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.bam
samtools index /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam
#./run_sv_caller.sh -A 1 -B 32 -z /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.discordant.bam -a /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam -b /m/cphg-RLscratch2/cphg_ratan/kk7t/target/sim4/sim100.all.sorted.bam -i ~/scratch/share/samtools-0.1.19/samtools -r /m/cphg-RLscratch/cphg-RLscratch/ar7jq/StructuralVariation/NA12878/alignments/na12878.ns.splitters.bam -C ../data/bed/ceph18.b37.exclude.2014-01-15.bed.webarchive -D ../results/text/ignoreCHR.txt
#time ($SAMTOOLS view -F 3586 -b -o ../data/bams/disc.bam ../data/bams/aln_test.bam)
#time ($SAMTOOLS view -F 3586 -b -o ../data/bams/disc.bam $REMOTE_FILE) #34m
#time ($SAMTOOLS view -f 64 -b -o ../data/bams/aln1s_100.bam ../data/bams/disc.bam) #4m
#time ($SAMTOOLS view -f 128 -b -o ../data/bams/aln2s_100.bam ../data/bams/disc.bam) # 4m

# run set cover
#./runSetCover7_indiv_temp.sh
#python formatBP_chrNum_v2.py ${TARGETNAME}.bedpe
#./run_classify.sh 25
#./run_PE_RD.sh $7 $8
