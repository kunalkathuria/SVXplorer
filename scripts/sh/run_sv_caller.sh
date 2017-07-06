#!/bin/bash
set -eux

READTHRESH=20
MATCHRATIO=1
NMATCHRATIO=0.95
NMATCH_PCT=0
CALC_THRESH=500000
RDS=0
VAR_RATE=250
SLOP=0
BP_MARGIN=20
SS=10
SR_MARGIN=250
RO=0.7
REF_RATE=5
DTHRESH=4
S_RISK=0
SIG_MULT=5
SIG_BOUND=2
MIN_CS=4
USC=1
SAM=samtools
BAM=
BAM_NS=
S_PATH=
RD_FILE="RDSegments.txt"

usage()
{
cat << EOF
usage: $0 [options]

Options (please use short options for now):
    -h 				(this menu)
    -a|--bamfile 		(*REQUIRED: position-sorted bam file)
    -b|--bam_ns 		(*REQUIRED: name-sorted bam file)
    -c|--read_thresh 		(max number of discordant read mappings a single PE fragment mate can have for it not to be ignored in uniqeness-based set cover, default = 20)
    -d|--match_thresh 		(min threshold of ratio of a given secondary alignment's score with AS of primary alignment needed for this secondary alignment to be considered/used in set cover, def = 1)
    -e|--n_match_thresh 	(same as above for number of base pair matches in alignment, instead of AS, def = .95)
    -f|--n_pct_thresh 		(min ratio of matching bases in a given alignment to mean read length required for read to be used in analysis, default =  0)
    -g|--calc_thresh 		(max entries of bamfile used to ascertain stats like coverage, mean insert length etc., def = 500,000)
    -i|--sam_path 		(complete path to samtools command, including /samtools, if not on path)
    -j|--bp_margin 		(in deciding breakpoint start, stop this is the minimum margin added to either side of breakpoint to account for sequencing errors etc. unless the calculated margin around breakpoints is anyway larger, def = 20 --> 10 added to either side)
    -k|--splitter_rate		(to save time, only 1 split read mapping within these many bases of others is used to improve breakpoint precision, since only 1 is really required if coming from same junction/breakpoint, def = 250)
    -l|--slop 			(in SV determination using PE mappings, this extra slop is used in calculating overlap as necessary for whatever reason, not needed ordinarily, def = 0)
    -m|--rd_signal 		(SET to 1 if using -x; whether or not to use read depth signal, e.g. if CN file is not available set this to 0, def = 0)
    -n|--rec_overlap		(min reciprocal overlap needed with CN region to boost/give more weight to intermediate PE-claimed SV)
    -o|--var_rate		(2 clusters are compared for overlap within this much max distance between breakpoint start locations, while consolidating clusters into variants; the lower this margin is, the faster the code will be; value should be around high-percentile of breakpoint region width and so depends upon coverage, def = 250)
    -p|--refresh_rate	 	(too specific, broadly -- controls how often variant list cleaned up. See ClassifyVariants.py, def = 5)
    -q|--disj_thresh		(number of disjoint/unique fragment mappings required in variant set to pick set as claimed SV, def = 4)
    -r|--splitter_path		(*REQUIRED: path to split-reads bam file)
    -s|--splitter_slop		(extra slop in determining SR overlap with existing PE variants, shouldn't be needed ordinarily, def = 0)
    -t|--splitter_risk	(too specific, broadly if this is set may gain both TPs and FPs, def = 0. See add_SR.py)
    -u|--min_cluster		(the minimum number of reads supporting a cluster required for it to be considered in the analysis, def = 5)
    -v|--sig_bound		(in forming breakpoint region, the value of mean_insert_length + sig_bound*sig_insert_length - 2*mean_RDL - (outer dist b/w min and max location of cluster read) is added as a margin to the "exact" breakpoint given by end-point of last cluster read in direction of orientation. This is usually 3 for Poisson etc., def = 2)
    -w|--sig_mult		(in forming clusters, a read's alignment location needs to fall within mean_IL + sig_mult*sigma_IL - 2*RDL of given cluster to be part of it as a necessary condition, def = 5)
    -x				(SET -m to 1 if using this; path to file containing CNVs)
    -y				(binary switch; set to 0 if wish to see only those variants that are supported uniquely by fragments, def = 1)	
EOF
}

while getopts “ha:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:” OPTION
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
	    ?)
	    echo "Option not recognized."
	    exit
	    ;;
     esac
done

if [[ -z $BAM ]] || [[ -z $BAM_NS ]] || [[ -z $SAM ]] || [[ -z $S_PATH ]]
then
     usage
     echo "Please provide paths to input BAM files as well as SAMTOOLS: required option fields are -a, -b, -r."
     exit 1
fi

time ($SAM view -F 3586 -b -o ../data/bams/discordants.bam $BAM_NS) #34m
time ($SAM view -f 64 -b -o ../data/bams/aln1s.bam ../data/bams/discordants.bam) #4m
time ($SAM view -f 128 -b -o ../data/bams/aln2s.bam ../data/bams/discordants.bam) # 4m

./form_clusters.sh $BAM_NS $BAM $READTHRESH $MATCHRATIO $NMATCHRATIO $CALC_THRESH $NMATCH_PCT $MIN_CS $SIG_MULT $BP_MARGIN $SIG_BOUND $SAM
./run_PE.sh $MIN_CS $VAR_RATE $SLOP $REF_RATE $DTHRESH $USC $BAM # can replace $MIN_CS here with another value and run these 2 steps with new threshold rather than rerun whole cluster formation sequence. In this case, simply comment out previous line and run sv_caller.
./run_RD_SR.sh $S_PATH $S_RISK $SS $SR_MARGIN $RO $DTHRESH $RDS $RD_FILE $USC
