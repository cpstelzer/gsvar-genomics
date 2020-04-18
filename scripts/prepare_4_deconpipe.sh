#!/bin/bash

# Delete files from previous decontamination_pipeline runs
# before a new run

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

K2OUT="/aqtrim/unmapped/kraken2/"
ROTIREADS="/rotireads/"

cd ..  # escape "contaminants" folder

while read CLONE SAMPLENO FRAGSIZE

do

    echo "Currently processing  $CLONE"	

    K2FOLDER="$CLONE$K2OUT"
    R_FOLDER="$CLONE$ROTIREADS"

    rm -r $K2FOLDER
    rm -r $R_FOLDER
   
done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
