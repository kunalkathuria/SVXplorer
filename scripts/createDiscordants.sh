#!/usr/bin/env bash
# run as ./createDiscordants.sh BAMFILE OUTPUT_DISC_FILE_PATH
BAM=$1
OUTPATH=$2
SAMTOOLS=samtools
THREADS=$3

$SAMTOOLS view -@$THREADS -F 3842 -b -o $OUTPATH $BAM
