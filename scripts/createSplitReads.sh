#!/usr/bin/env bash
# run as ./createSplitReads.sh ALMTFILE(SORTED BY POS OR NAME) OUTPUT_FILE_PATH PATH_TO_LUMPY N_THREADS
BAM=$1
OUTFILE=$2
LUMPY=$3
SAMTOOLS=samtools
THREADS=$4

$SAMTOOLS view -@$THREADS -h $BAM \
        | $LUMPY/scripts/extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - \
        > $OUTFILE
