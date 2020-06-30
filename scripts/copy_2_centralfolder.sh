#!/bin/bash

# Collecting clone-specific data and store it to a central folder

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="/finalreads/"
OUTMAIN="results/finalreads/fq/"  # target location

cd ..  # escape "SCRIPTS" folder
# mkdir $OUTMAIN


while read CLONE SAMPLENO FRAGSIZE

do
  	
    echo "Currently processing  $CLONE"	
    # OHJ13_mbl2016_to_VBCFpol.sorted.bam
    INFILE1="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R1.fq.gz""
    INFILE2="$CLONE$INMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R2.fq.gz""
    OUTFILE1="$OUTMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R1.fq.gz""
    OUTFILE2="$OUTMAIN$CLONE"_"$SAMPLENO".roti-mito.dedup.R2.fq.gz""
    cp $INFILE1 $OUTFILE1
    cp $INFILE2 $OUTFILE2


done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS


