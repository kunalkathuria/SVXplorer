#!/usr/bin/env bash
# run as ./createSplitReads.sh BAMFILE OUTPUT_DISC_FILE_PATH $PATH_TO_LUMPY
BAM=$1
OUTPATH=$2
LUMPY=$3
SAMTOOLS=samtools
THREADS=$4

$SAMTOOLS view -@$THREADS -h $BAM \
        | $LUMPY/scripts/extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - \
        > $OUTPATH
