#!/bin/bash

# Collecting clone-specific data and store it to a central folder

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="/finalmapping/"
OUTMAIN="results/finalmapping/"  # target location

cd ..  # escape "SCRIPTS" folder
mkdir $OUTMAIN


while read CLONE SAMPLENO FRAGSIZE

do
  	
    echo "Currently processing  $CLONE"	
    # OHJ13_mbl2016_to_VBCFpol.sorted.bam
    INFILE="$CLONE$INMAIN$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam""
    OUTFILE="$OUTMAIN$CLONE"_"$SAMPLENO"_to_VBCFpol.sorted.bam""
    cp $INFILE $OUTFILE


done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS


