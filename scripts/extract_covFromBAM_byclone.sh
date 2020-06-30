#!/bin/bash

# Extracting mean coverage from BAM file after combining all alignments of the same clone (see Notebook p. 113, 3.6.2020)

OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

INMAIN="finalmapping"
ASSEMBLY="VBCFpol"

cd ..

# Make a list with unique clone names
while read CLONE SAMPLENO FRAGSIZE
do
    echo "$CLONE" >> clones1.tmp
done < SCRIPTS/datasets/OHJ_all.csv

cat clones1.tmp | sort | uniq > clones2.tmp

while read CLONE
do
    
    echo "Currently processing  $CLONE"
    
    INFOLDER="$CLONE"/"$INMAIN"/
    INBAMS="*$ASSEMBLY.sorted.bam"
 
    cd $INFOLDER
    samtools merge ../../temp.merged.bam $INBAMS
    cd ../..

    AVGCOV=$( samtools depth temp.merged.bam | awk '{sum+=$3} END {print sum/NR}')

    COL[1]=$CLONE
    COL[2]=$AVGCOV

    IFS=$OLDIFS
     echo ${COL[*]} >> SCRIPTS/cov_from_bam.report.byclone.csv
    IFS=$NEWIFS

    rm temp.merged.bam

done < clones2.tmp

rm clones1.tmp clones2.tmp

IFS=$OLDIFS

