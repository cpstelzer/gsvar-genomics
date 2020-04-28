#!/bin/bash

# Script for collecting information on bbduk-trimming results
# See notebook (p.34, 5.3.2020)



OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="/aqtrim/"

cd ..

COL[1]="Collected Information from bbduk-trimming reports"
echo ${COL[*]} >> SCRIPTS/aqtrim.report.csv
COL[1]="Library"
COL[2]="total_reads"
COL[3]="recovered_reads"
COL[4]="total_bases"
COL[5]="recovered_bases"
echo ${COL[*]} >> SCRIPTS/aqtrim.report.csv


while read CLONE LIBRARY FRAGSIZE

do
     INFILE="$CLONE$INMAIN$CLONE"_"$LIBRARY.report.txt"
     TEMPFILE="$CLONE$INMAINgrepped.tmp"
     
     TOTALREADS=$(grep -Po '(?<=Input:\s{18}\t)[0-9]*(?=\sreads)' $INFILE)
     RECOVEREDREADS=$(grep -Po '(?<=Result:\s{17}\t)[0-9]*(?=\sreads)' $INFILE)
     TOTALBASES=$(grep -Po '[0-9]*(?=\sbases\.)' $INFILE)
     cat $INFILE | grep "Result:" > $TEMPFILE
     RECOVEREDBASES=$(grep -Po '[0-9]*(?=\sbases)' $TEMPFILE) 
     rm $TEMPFILE
     
     COL[1]="$CLONE"_"$LIBRARY"
     COL[2]=$TOTALREADS
     COL[3]=$RECOVEREDREADS
     COL[4]=$TOTALBASES
     COL[5]=$RECOVEREDBASES
     echo ${COL[*]} 
     
     IFS=$OLDIFS
     echo ${COL[*]} >> SCRIPTS/aqtrim.report.csv
     IFS=$NEWIFS

done < SCRIPTS/datasets/OHJ_all.csv

IFS=$OLDIFS
