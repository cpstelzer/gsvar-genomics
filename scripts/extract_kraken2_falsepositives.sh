#!/bin/bash

# Pipeline zur Extraktion der "false_positives" der kraken2 Analyse
# Die false positive rate entspricht der % der ans Referenzgenom mappenden Reads,
# die ebenso ein positives Signal bei der Kraken2 Analyse ergibt
# Note: Saving CSV file does not work properly --> save screen output instead

# CP Stelzer 19.5.2020



cd ..  # escape "contaminants" folder

KRAKENPATH="contaminants/kraken2/false_positives3/"
OUTREPORT="contaminants/kraken2.false_positives.csv"

COL[1]="clone"
COL[2]="library"    
COL[3]="std_db"
COL[4]="proto_db"
echo ${COL[*]} >> $OUTREPORT

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

while read CLONE SAMPLENO FRAGSIZE

do # loop through clones/libraries
    
    FILE_STD="$KRAKENPATH$CLONE"_"$SAMPLENO".mapped.std_db.log"" # standard database (-human)
    FILE_PROTO="$KRAKENPATH$CLONE"_"$SAMPLENO".mapped.proto_db.log"" # protozoan database
    FPR_STD=$(sed -n '2s/ *//p' $FILE_STD | cut -f 1)    
    FPR_PROTO=$(sed -n '2s/ *//p' $FILE_PROTO | cut -f 1)
    echo "$CLONE $SAMPLENO $FPR_STD $FPR_PROTO"

    IFS=$OLDIFS    
    COL[1]=$CLONE
    COL[2]=$SAMPLENO    
    COL[3]=$FPR_STD
    COL[4]=$FPR_PROTO
    echo ${COL[*]} >> "contaminants/deconreport.csv"
    IFS=$NEWIFS 


done < SCRIPTS/datasets/OHJ_all.csv

IFS=$OLDIFS
