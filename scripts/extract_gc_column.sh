#!/bin/bash

# Extract column in which GC values are stored 
# (after running 'read_gc_content.sh')

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

OUTDIR="gc_dist/"
mkdir -p $OUTDIR

while read CLONE SAMPLENO FRAGSIZE

    do
    echo "Currently processing $CLONE $SAMPLENO"
    INFILE="$OUTDIR$CLONE"_"$SAMPLENO".roti-mito.dedup.gc.tab""
    OUTFILE="$OUTDIR$CLONE"_"$SAMPLENO".gc.csv""
    awk '{print $2}' $INFILE > $OUTFILE

done < OHJ_all.csv


IFS=$OLDIFS
