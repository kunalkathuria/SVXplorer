set -eux

BWA=~/scratch/share/bwa-0.7.12/bwa
SAMTOOLS=/m/cphg-RLscratch/cphg-RLscratch/share/samtools-0.1.19/samtools
WGSIM=~/scratch/share/wgsim/wgsim
SVSIM=./sv_sim/SVsim
SVCommands=./sv_sim/SVcommands_test_1SV.sim
REFERENCE=/m/cphg-RLscratch/cphg-RLscratch/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta
TARGET=~/other/target/genome_rearranged.fasta
READ1=../data/fqs/read_1.fq
READ2=../data/fqs/read_2.fq
REMOTE_FILE=~/homozygous/alignments/id.ns.alignment.bam

MEAN=500
STDEV=50
RDL=100
NBP=3100000000
COVERAGE=50
NREADS=$(( (COVERAGE * NBP) / (2 * RDL) ))
NSIMS=1

# IF TARGET IS NOT PREPARED
# prepare simualted target
#python $SVSIM -i $SVCommands -r $REFERENCE -o $TARGETNAME -W -d -n $NSIMS

# prep fake reads from the target
#$WGSIM -N$NREADS -1 $RDL -2 $RDL -d $MEAN -s $STDEV $TARGET $READ1 $READ2

# bwa alignment
#$BWA index $REFERENCE
#$BWA mem -R '@RG\tID:foo\tSM:bar' -a -Y -t 32 $REFERENCE $READ1 $READ2 \
#| $SAMTOOLS view -S -b - \
#> ../data/bams/test/test_target_pe.bam

#$SAMTOOLS index ../data/bams/test/test_target_pe.bam
#$SAMTOOLS sort -@ 32 -n ../data/bams/test/test_target_pe.bam ../data/bams/test/test_target_pe.ns

/m/cphg-RLstore/cphg-RLstore/kk7t/sv_caller/other_tools/lumpy/code/lumpy-sv-0.2.11/run_lumpy.sh

./run_sv_caller.sh -a ../data/bams/test/test_target_pe.bam -b ../data/bams/test/test_target_pe.ns.bam -i $SAMTOOLS -r ../data/bams/test/test_target_splitters.ns.bam
#time ($SAMTOOLS view -F 3586 -b -o ../data/bams/discordants.bam $REMOTE_FILE) #34m
#time ($SAMTOOLS view -f 64 -b -o ../data/bams/aln1s.bam ../data/bams/discordants.bam) #4m
#time ($SAMTOOLS view -f 128 -b -o ../data/bams/aln2s.bam ../data/bams/discordants.bam) # 4m

# run set cover
#./form_clusters_d.sh
#./run_PE_d.sh 4
#./run_RD_SR_d.sh
