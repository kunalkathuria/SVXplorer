#!/usr/bin/env bash
# run as ./createSplitReads.sh BAMFILE OUTPUT_DIR OUTPUT_FILE_NAME PATH_TO_LUMPY N_THREADS
BAM=$1
OUTDIR=$2
OUTFILE=$3
LUMPY=$4
SAMTOOLS=samtools
THREADS=$5

$SAMTOOLS view -@$THREADS -h $BAM \
        | $LUMPY/scripts/extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - \
        > $OUTDIR/temp.bam
$SAMTOOLS view -@$THREADS -q 1 -b -h -o $OUTFILE $OUTDIR/temp.bam
rm $OUTDIR/temp.bam
