#!/bin/bash
set -eux

READTHRESH=20
ABS_AS=0
MAP_THRESH=1
MAP_THRESH_U=-1
MQ_SR=20
MATCHRATIO=.95
NMATCHRATIO=0
NMATCH_PCT=0
CALC_THRESH=1000000
SIG_THRESH=.9985
RDS=0
VAR_RATE=450
SLOP=0
BP_MARGIN=30
SS=10
SR_MARGIN=0
RO=0.7
REF_RATE=5
DTHRESH=3
S_RISK=1
SIG_MULT=3 #obsolete
SIG_BOUND=3
MIN_CS=3
MIN_CS_SR=$DTHRESH
USC=1
SAM=samtools
BAM=
BAM_NS=
S_PATH=
RD_FILE="../results/text/cnvs.err.txt"
DISCORDANTS=
THREADS=1
DISC_IS_NS=1
IGNORE_BED="none"
IGNORE_CHR="none"
CLBREAK_QUANTILE=.9999
LIB_MULT=2.0
DISC_ENHANCER3=2.0
AS_STAT_QUANTILE=.999
DISC_ENHANCER=1.67
CL_BRK=1
LIB_RDL_GAP=1.0
RD_VARNUM=3000000
RD_FRAGNUM=100000000
RD_MARGIN=0
SR_INS_SEPARATION=30
RD_SLOP=250
RD_MQU=10
MINSR_PE_SUPP=2
DEL_CN_SUPP_THRESH=.75 
DUP_CN_SUPP_THRESH=1.25

usage()
{
cat << EOF
usage: $0 [options]

Options (please use short options for now):
    -h 				(this menu)
    -a|--bamfile 		(*REQUIRED: path of position-sorted bam file)
    -b|--bam_ns 		(*REQUIRED: path of name-sorted bam file)
    -z|--disc			(*REQUIRED: path of BAM file containing all discordant pairs extracted from above. See README) 
    -A|ns			(Set this to 0 if discordant file above is not name-sorted; def = 1)
    -B|--threads		(Number of available cores; def = 1)
    -C|--ignoreB		(BED file in format "chr start stop", w/o ANY headers, containing reference regions to ignore in forming variants)
    -D|--ignoreChr		(file containing names of chromosomes, 1 per line,  whose alignments should be completely ignored)
    -c|--read_thresh 		(max number of discordant read mappings a single PE fragment mate can have for it not to be ignored in uniqeness-based set cover, default = 20)
    -d|--match_thresh 		(min threshold of ratio of a given secondary alignment's score with AS of primary alignment needed for this secondary alignment to be considered/used in set cover, def = .95)
    -e|--n_match_thresh 	(same as above for number of base pair matches in alignment, instead of AS, def = 0)
    -f|--n_pct_thresh 		(min ratio of matching bases in a given alignment to read length required for read to be used in analysis, default =  0)
    -g|--calc_thresh 		(max entries of bamfile used to ascertain stats like coverage, mean insert length etc., def = 1,000,000)
    -i|--sam_path 		(complete path to "samtools" command, including /samtools, if not on path)
    -m|--rd_signal 		(SET to 1 if using -x; whether or not to use read depth signal, e.g. if CN file is not available set this to 0, def = 0)
    -n|--rec_overlap		(min reciprocal overlap needed with CN region to boost/give more weight to intermediate PE-supported SV)
    -q|--disj_thresh		(number of unique fragment mappings required in variant set to pick set as claimed SV, def = 4)
    -r|--splitter_path		(*REQUIRED: path to split-reads bam file)
    -s|--splitter_slop		(margin used in differentiating split-read break points from one another; def = 10)
    -u|--min_cluster		(the minimum number of reads supporting a cluster required for it to be considered in the analysis, def = 3)
    -x				(SET -m to 1 if using this; path to file containing CNVs)
    -y				(binary switch; set to 0 if wish to see only those variants whose supporting fragments do not map at all to other variants, def = 1)	
EOF
}

while getopts “ha:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:” OPTION
do

     case $OPTION in
	    h)
	    usage
	    exit 1
	    ;;	
    	    a)
	    BAM="$OPTARG"
	    ;;
	    b)
	    BAM_NS="$OPTARG"
	    ;;
	    c)
	    READTHRESH="$OPTARG"
	    ;;
	    d)
	    MATCHRATIO="$OPTARG"
	    ;;
	    e)
	    NMATCHRATIO="$OPTARG"
	    ;;
	    f)
	    NMATCH_PCT="$OPTARG"
	    ;;
	    g)
	    CALC_THRESH="$OPTARG"
	    ;;
	    i)
            SAM="$OPTARG"
            ;;
	    j)
	    BP_MARGIN="$OPTARG"
	    ;;
	    k)
	    SR_MARGIN="$OPTARG"
	    ;;
	    l)
	    SLOP="$OPTARG"
	    ;;
	    m)
            RDS="$OPTARG"
            ;;
	    n)
	    RO="$OPTARG"
	    ;;
	    o)
	    VAR_RATE="$OPTARG"
	    ;;
	    p)
	    REF_RATE="$OPTARG"
	    ;;
	    q)
	    DTHRESH="$OPTARG"
	    ;;
	    r)
	    S_PATH="$OPTARG"
	    ;;
	    s)
	    SS="$OPTARG"
	    ;;
	    t)
	    S_RISK="$OPTARG"
	    ;;
	    u)
	    MIN_CS="$OPTARG"
	    ;;
	    v)
	    SIG_BOUND="$OPTARG"
	    ;;
	    w)
	    SIG_MULT="$OPTARG"
	    ;;
	    x)
	    RD_FILE="$OPTARG"
	    ;;
	    y)
	    USC="$OPTARG"
	    ;;
	    z)
            DISCORDANTS="$OPTARG"
            ;;
	    A)
            DISC_IS_NS="$OPTARG"
            ;;
	    B)
            THREADS="$OPTARG"		
	    ;;
	    C)
	    IGNORE_BED="$OPTARG"
	    ;;
	    D)
	    IGNORE_CHR="$OPTARG"
	    ;;
	    E)
            CL_BRK="$OPTARG"
            ;;
	    ?)
	    echo "Option not recognized."
	    exit
	    ;;
     esac
done

if [[ -z $BAM ]] || [[ -z $BAM_NS ]] || [[ -z $S_PATH ]]
then
     usage
     echo "Please provide paths to input BAM files as well as SAMTOOLS: required option fields are -a, -b, -z, -r."
     exit 1
fi

#IF disc BAM not ready
#time ($SAM view -@32 -F 3586 -b -o /m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR505885/disc.all.bam $BAM_NS) #34m
#DISCORDANTS=/m/cphg-RLscratch2/cphg_ratan/kk7t/target/NA12878/SRR3397076/disc.all.bam

#Remove next 2 comments when done testing
time ($SAM view -@32 -f 64 -b -o ../data/bams/aln1s_u.bam $DISCORDANTS) #4m
time ($SAM view -@32 -f 128 -b -o ../data/bams/aln2s_u.bam $DISCORDANTS) # 4m
 
if [ $DISC_IS_NS -eq 1 ]
then
	mv ../data/bams/aln1s_u.bam ../data/bams/aln1s.bam # change back to al1ns_u.bam etc
	mv ../data/bams/aln2s_u.bam ../data/bams/aln2s.bam
else
	$SAM sort -n -@ $THREADS ../data/bams/aln1s_u.bam ../data/bams/aln1s
	$SAM sort -n -@ $THREADS ../data/bams/aln2s_u.bam ../data/bams/aln2s
fi

./form_clusters.sh $BAM_NS $BAM $READTHRESH $MATCHRATIO $NMATCHRATIO $CALC_THRESH $NMATCH_PCT $MIN_CS $SIG_MULT $BP_MARGIN $SIG_BOUND $SAM $IGNORE_BED $IGNORE_CHR $MAP_THRESH $ABS_AS $MAP_THRESH_U $SIG_THRESH $AS_STAT_QUANTILE $CLBREAK_QUANTILE $LIB_MULT $DISC_ENHANCER3 $DISC_ENHANCER $CL_BRK $LIB_RDL_GAP
./run_PE.sh $MIN_CS $VAR_RATE $SLOP $REF_RATE $DTHRESH $USC $BAM $MAP_THRESH_U $SIG_MULT $RD_VARNUM $RD_FRAGNUM $DEL_CN_SUPP_THRESH $DUP_CN_SUPP_THRESH # can replace $MIN_CS here with another value and run these 2 steps with new threshold rather than rerun whole cluster formation sequence. In this case, simply comment out previous line and run sv_caller.
./run_RD_SR.sh $S_PATH $S_RISK $SS $SR_MARGIN $RO $DTHRESH $RDS $RD_FILE $USC $IGNORE_CHR $MAP_THRESH_U $MQ_SR $MIN_CS_SR $BAM $RD_VARNUM $RD_FRAGNUM $RD_MARGIN $SR_INS_SEPARATION $RD_SLOP $RD_MQU $MINSR_PE_SUPP $DEL_CN_SUPP_THRESH $DUP_CN_SUPP_THRESH
