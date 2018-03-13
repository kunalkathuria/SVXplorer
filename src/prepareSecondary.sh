#!/usr/bin/env bash
THREADS=32
BAMFILE1=$1
BAMFILE2=$2
WORK_DIR=$3

samtools view -@ $THREADS -F 1294 -b -o $WORK_DIR/disc1.bam $BAMFILE1
samtools view -@ $THREADS -F 1294 -b -o $WORK_DIR/disc2.bam $BAMFILE2
samtools merge -nf $WORK_DIR/discordants.bam $WORK_DIR/disc1.bam $WORK_DIR/disc2.bam 
