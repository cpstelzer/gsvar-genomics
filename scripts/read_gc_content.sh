##!/bin/bash

# Script to extract the gc content of each read from fastq-files
# Run on leo4 in the same folder as the data
# This file has to be activated in the bioconda environment
# called from: submit_read_gc.sge

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

OUTDIR="gc_dist/"
mkdir -p $OUTDIR

while read CLONE SAMPLENO FRAGSIZE

    do
    INFILE1="$CLONE"_"$SAMPLENO".roti-mito.dedup.R1.fq.gz""
    seqtk sample -s100 $INFILE1 10000 > temp.IN1.fq
    INFILE2="$CLONE"_"$SAMPLENO".roti-mito.dedup.R2.fq.gz""
    seqtk sample -s100 $INFILE2 10000 > temp.IN2.fq
    OUTFILE="$OUTDIR$CLONE"_"$SAMPLENO".roti-mito.dedup.gc.tab""

    seqkit fx2tab --name --only-id --gc temp.IN1.fq > $OUTFILE
    seqkit fx2tab --name --only-id --gc temp.IN2.fq >> $OUTFILE

done < OHJ_all.csv


IFS=$OLDIFS

