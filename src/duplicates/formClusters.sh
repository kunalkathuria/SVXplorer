#!/bin/bash
#set -eux

if [ -d ./results ]
then
        rm -r ./results
fi
mkdir ./results

BAMFILE=
IGNORE_POS="none"
IGNORE_CHR="none"
MAP_THRESH=
AS_THRESH=0
MIN_CLUSTER_SIZE=
SAM=samtools
DISC_IS_NS=0
THREADS=32

usage()
{
cat << EOF
usage: $0 [options]

Options (please use short options for now):
    -h                          (this menu)
    -b                          (*REQUIRED: path of position-sorted or name-sorted bam file)
    -m                          (*REQUIRED: Min number of alignments needed for cluster)
    -i                          (path to BED file to ignore alignments falling in these regions)
    -c                          (path to file containing list of chromosome/contig names, 1 per line, to be ignored)
    -t                          (*REQUIRED: mapping quality threshold)
    -s                          (alignment score threshold)

All results and intermediatte files stored in ./results

EOF
}

while getopts “hb:m:i:c:t:s:” OPTION
do

     case $OPTION in
            h)
            usage
            exit 1
            ;;
            b)
            BAMFILE="$OPTARG"
            ;;
            m)
            MIN_CLUSTER_SIZE="$OPTARG"
            ;;
            i)
            IGNORE_POS="$OPTARG"
            ;;
            c)
            IGNORE_CHR="$OPTARG"
            ;;
            t)
            MAP_THRESH="$OPTARG"
	    ;;
	    s)
	    AS_THRESH="$OPTARG"			
	    ;;
	    ?)
            echo "Option not recognized."
            exit
            ;;
     esac
done

if [[ -z $BAMFILE ]] || [[ -z $MIN_CLUSTER_SIZE ]] || [[ -z $MAP_THRESH ]]
then
     usage
     echo "Please provide paths to input BAM files as well as SAMTOOLS: required option fields are -a, -b, -z, -r."
     exit 1
fi

time ($SAM view -@32 -F 3586 -b -o ./results/disc.bam $BAMFILE) 
DISCORDANTS=./results/disc.bam

#Remove next 2 comments when done testing
time ($SAM view -@32 -f 64 -b -o ./results/aln1s_u.bam $DISCORDANTS) #4m
time ($SAM view -@32 -f 128 -b -o ./results/aln2s_u.bam $DISCORDANTS) # 4m

if [ $DISC_IS_NS -eq 1 ]
then
       mv ./results/aln1s_u.bam ./results/aln1s.bam # change back to al1ns_u.bam etc
       mv ./results/aln2s_u.bam ./results/aln2s.bam
else
       $SAM sort -n -@ $THREADS ./results/aln1s_u.bam > ./results/aln1s.bam
       $SAM sort -n -@ $THREADS ./results/aln2s_u.bam > ./results/aln2s.bam
fi

time ( python ReadDiscordants.prim.graph.py $BAMFILE $IGNORE_POS $IGNORE_CHR $MAP_THRESH $AS_THRESH) # 9 minutes on 50 x homozygous; last fields are map_thresh,AS_THRESH temporary
echo Sorting
sort -T ./results/ -k 2,2 -k 3,3n ./results/All_Discords_P.txt > ./results/All_Discords_P_S.txt
echo DoneSorting

time (python FormClusters.samechr.py ./results/bam_stats.txt $MIN_CLUSTER_SIZE) #> ../results/text/fc_time.txt # 50m

