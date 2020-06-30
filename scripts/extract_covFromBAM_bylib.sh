#!/bin/bash

# Extracting mean coverage from BAM file (see Notebook p. 113, 3.6.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS


INMAIN="finalmapping"
ASSEMBLY="VBCFpol"

cd ..


while read CLONE SAMPLENO ADDCOL

do
    INFILE="$CLONE/$INMAIN/$CLONE"_"$SAMPLENO"_to_"$ASSEMBLY".sorted.bam""
    AVGCOV=$( samtools depth $INFILE | awk '{sum+=$3} END {print sum/NR}')
    COL[1]=$CLONE
    COL[2]=$SAMPLENO
    COL[3]=$AVGCOV

    IFS=$OLDIFS
     echo ${COL[*]} >> SCRIPTS/cov_from_bam.report.bylib.csv
    IFS=$NEWIFS

done < SCRIPTS/datasets/OHJ_all.csv

IFS=$OLDIFS

