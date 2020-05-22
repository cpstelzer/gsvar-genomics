#!/bin/bash

# Counting reads and basepairs at intermediate steps of the preprocessing pipeline

# This script analyses read files (fastq) by:
# 1) unzipping fq.gz
# 2) convert from fastq to fasta
# 3) count #reads in R1 and R2
# 4) count #bases in R1 and R2
# 5) writing the results to a csv-file
# 6) zipping fq

# E.g. trimmed reads (after adaptor and quality trimming
#      rotireads (after initial contaminant removal)

# Secondary input file (preprocessing_steps.csv)
# rawdata;
# aqtrim;.trimmed
# rotireads;.rotireads
# finalreads;.roti-mito.dedup


OLDIFS=$IFS     # save the existing field separator
NEWIFS=";"
IFS=$NEWIFS

cd ..  # escape SCRIPTS folder


while read CLONE SAMPLENO FRAGSIZE
do
    echo "Currently processing  $CLONE"_" $SAMPLENO"

    NCOUNT=0
    COL[1]=$CLONE
    COL[2]=$SAMPLENO

    while read INPATH EXTENSION
    do

        echo "Extracting data from: $INPATH"

        NCOUNT=$NCOUNT+4

        READS1_GZ="$CLONE/$INPATH/$CLONE"_"$SAMPLENO$EXTENSION".R1.fq.gz""
        READS2_GZ="$CLONE/$INPATH/$CLONE"_"$SAMPLENO$EXTENSION".R2.fq.gz""
        READS1="$CLONE/$INPATH/$CLONE"_"$SAMPLENO$EXTENSION".R1.fq""
        READS2="$CLONE/$INPATH/$CLONE"_"$SAMPLENO$EXTENSION".R2.fq""

        gunzip $READS1_GZ $READS2_GZ
        
        FASTA_R1="temp.R1.fasta"
        FASTA_R2="temp.R2.fasta"
        sed -n '1~4s/^@/>/p;2~4p' $READS1 > $FASTA_R1  # Convert to FASTA-format
        sed -n '1~4s/^@/>/p;2~4p' $READS2 > $FASTA_R2
        N_READS_R1=$(grep -c '^>' $FASTA_R1)  # number of reads in multi-fasta
        N_READS_R2=$(grep -c '^>' $FASTA_R2)
        #(( COL[3] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
        (( COL[$NCOUNT-1] = $N_READS_R1 + $N_READS_R2 ))   # sum-up results of both mates
        TOTAL_BP_R1=$(grep -v '>' $FASTA_R1 | tr -d '\n' | wc -c) # total number of bases in fasta
        TOTAL_BP_R2=$(grep -v '>' $FASTA_R2 | tr -d '\n' | wc -c)
        #(( COL[4] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))
        (( COL[$NCOUNT] = $TOTAL_BP_R1 + $TOTAL_BP_R2 ))

        gzip $READS1 $READS2


    done < SCRIPTS/datasets/preprocessing_steps.csv

    echo ${COL[*]} >> "SCRIPTS/quantify.reads_bp.report.csv"

done < SCRIPTS/datasets/OHJ_all.csv


IFS=$OLDIFS
